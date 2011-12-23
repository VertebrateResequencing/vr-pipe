use VRPipe::Base;

class VRPipe::Manager extends VRPipe::Persistent {
    use Parallel::ForkManager;
    use Sys::CPU;
    use POSIX qw(ceil);
    
    our $DEFAULT_MAX_PROCESSES = Sys::CPU::cpu_count();
    
    has '_running' => (is => 'rw',
                       isa => 'Bool',
                       traits => ['VRPipe::Persistent::Attributes'],
                       default => 0);
    
    has '_run_start' => (is => 'rw',
                         isa => Datetime,
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    has '_instance_run_start' => (is => 'rw',
                                  isa => Datetime);
    
    has 'global_limit' => (is => 'rw',
                           isa => PositiveInt);
    
    has '_created_pipelines' => (is => 'rw',
                                 isa => 'HashRef');
    
    # public getters for our private attributes
    method running {
        return $self->_running;
    }
    method run_start {
        return $self->_run_start;
    }
    
    __PACKAGE__->make_persistent();
    
    around get (ClassName|Object $self: Persistent :$id?, PositiveInt :$global_limit = 500) {
        # our db-storage needs are class-wide, so we only have 1 row in our
        # table
        return $self->$orig(id => 1, global_limit => $global_limit);
    }
    
    method setups (Str :$pipeline_name?) {
        my $rs = $self->result_source->schema->resultset('PipelineSetup');
        my @setups;
        while (my $ps = $rs->next) {
            my $p = $ps->pipeline;
            if ($pipeline_name) {
                next unless $p->name eq $pipeline_name;
            }
            
            # before spooling, whilst we are still in a single process, make
            # sure that all pipelines have had their step members created using
            # the steps method. *** steps() currently gives a memory leak, so
            # that's the other reason we are careful to only call it once per
            # pipeline
            my $created = $self->_created_pipelines;
            unless (exists $created->{$p->id}) {
                $p->steps;
                $created->{$p->id} = 1;
                $self->_created_pipelines($created);
            }
            
            push(@setups, $ps);
        }
        
        return @setups;
    }
    
    method trigger (ArrayRef[VRPipe::PipelineSetup] :$setups?, PositiveInt :$max_processes?) {
        $setups ||= [$self->setups];
        $max_processes ||= $DEFAULT_MAX_PROCESSES;
        
        # spool through multiple setups at once
        my $fm = Parallel::ForkManager->new($max_processes);
        
        # we'll collect data from our spools
        my %retrieved_responses = ();
        $fm->run_on_finish (
            sub {
                my $data_structure_reference = $_[5];
                if (defined($data_structure_reference)) {
                    $retrieved_responses{$$data_structure_reference} = 1;
                }
                else {
                    $self->throw("spool child failed to return anything");
                }
            }
        );
        
        foreach my $setup (@$setups) {
            if ($setup->active) {
                $fm->start and next; # fork
                
                my $all_done = $self->spool($setup);
                
                $fm->finish(0, \$all_done); # exit in the child process
            }
            else {
                #if (time() - $setup->last_updated > 7776000) {
                #    # it's been inactive for ~3months; perhaps we want to trash old stuff?
                #    $self->expunge($setup);
                #    # expunges all Persistent entries in the Peristent db that were created
                #    # during the course of this pipeline setup, that are not used by another
                #    # pipeline setup
                #}
            }
        }
        $fm->wait_all_children;
        
        my @results = keys %retrieved_responses;
        if (@results == 1 && $results[0] == 1) {
            return 1;
        }
        return 0;
    }
    
    method spool (VRPipe::PipelineSetup $setup) {
        # This first does VRPipe::Persistent::capture_start and when done will
        # call capture_stop, and will turn the array into a PersistentArray
        # which it will store in $set->children, merging the array with any
        # existing PersistentArray. This capture_stop&merge code will also run
        # during destruction of $man to protect against being killed before
        # run() completes normally. Under normal circumstances it will also
        # touch() every child of $pip, to prevent expunge_old destroying our
        # lives.
        # It also calls $self->running(1) which sets a boolean to true in the
        # persistent storage, along with the time we called this. Before
        # returning from run(), or on destruction, we set running(0) if the
        # timestamp was ours. If another process comes along and calls ->run
        # while we are running(), then the other process will block for some
        # reasonable amount of time until running is 0. If it takes too long,
        # then it will stop blocking and update the timestamp on running() and
        # continue with its own run().
        #
        # Next, it uses $set->data_source (having checked it is one of the
        # $set->module's valid ones) to loop through each data input. For each
        # data_source element it creates an instance of the $set->module. It
        # uses the ->steps method on that to set up
        # VRPipe::PipelineManager::Step objects and uses those to track state
        # as it does ->run_step($name) etc. Note that this process is only
        # creating Step and Submission objects; we are not calling ->submit
        # on the Submissions, so run() doesn't itself directly cause anything
        # to be executed on the farm.
        #
        # As it loops through everything it keeps track of various things like
        # number of data_source elements, steps completed and still to go per
        # element and overall, and number of jobs
        # completed/running/to-be-run/failed per step, per element and
        # overall. It returns the total number of jobs still to go (those
        # currently running, or scheduled to run later), so if this is 0 then
        # this pipeline is finished for now at least.
        
        # capture start
        
        # touch every child of setup
        
        my $setup_id = $setup->id;
        my $pipeline = $setup->pipeline;
        my @step_members = $pipeline->step_members;
        my $num_steps = scalar(@step_members);
        
        my $datasource = $setup->datasource;
        my $output_root = $setup->output_root;
        $self->make_path($output_root);
        my $all_done = 1;
        my $incomplete_elements = 0;
        my $limit = $self->global_limit;
        
        # We don't pass a limit to incomplete_element_states since that could
        # return all estates where all their submissions have failed 3 times.
        # To avoid doing wasted cycles of the estates loop we do make use of
        # limit though. The real limiting happens in handle_submissions()
        my $estates = $datasource->incomplete_element_states($setup);
        
        $self->debug("Spool for setup ".$setup->id." got ".scalar(@$estates)." incomplete element states, given a limit of $limit");
        
        foreach my $estate (@$estates) {
            my $element = $estate->dataelement;
            my $completed_steps = $estate->completed_steps;
            $self->debug("setup $setup_id | incomplete element loop $incomplete_elements for element id ".$element->id." which has completed $completed_steps steps of $num_steps");
            next if $completed_steps == $num_steps; # (rare, if ever?)
            
            $incomplete_elements++;
            
            my %previous_step_outputs;
            my $already_completed_steps = 0;
            foreach my $member (@step_members) {
                my $step_number = $member->step_number;
                my $state = VRPipe::StepState->get(stepmember => $member,
                                                   dataelement => $element,
                                                   pipelinesetup => $setup);
                
                my $step = $member->step(previous_step_outputs => \%previous_step_outputs, step_state => $state);
                if ($state->complete) {
                    $self->_complete_state($step, $state, $step_number, $pipeline, \%previous_step_outputs);
                    $already_completed_steps++;
                    
                    if ($already_completed_steps > $completed_steps) {
                        $estate->completed_steps($already_completed_steps);
                        $completed_steps = $already_completed_steps;
                        $estate->update;
                    }
                    
                    if ($already_completed_steps == $num_steps) {
                        # this element completed all steps in the pipeline
                        $incomplete_elements--;
                    }
                    
                    next;
                }
                
                # have we previously done the dispatch dance and are currently
                # waiting on submissions to complete?
                my @submissions = $state->submissions;
                if (@submissions) {
                    my $unfinished = $self->unfinished_submissions(submissions => \@submissions);
                    unless ($unfinished) {
                        my $ok = $step->post_process();
                        if ($ok) {
                            $self->_complete_state($step, $state, $step_number, $pipeline, \%previous_step_outputs);
                            next;
                        }
                        else {
                            # we warn instead of throw, because the step may
                            # have discovered its output files are missing and
                            # restarted itself
                            $self->warn("submissions completed, but post_process failed");
                        }
                    }
                    else {
                        # don't count this toward the limit if all the
                        # submissions have failed 3 times *** max_retries needs to be consistent and set between trigger and handle
                        my $fails = 0;
                        foreach my $sub (@submissions) {
                            next unless $sub->failed;
                            if ($sub->retries >= 3) {
                                $fails++;
                            }
                        }
                        
                        if ($fails == @submissions) {
                            $incomplete_elements--;
                        }
                    }
                }
                else {
                    # this is the first time we're looking at this step for
                    # this data member for this pipelinesetup
                    my $completed;
                    try {
                        $completed = $step->parse();
                    }
                    catch ($err) {
                        warn $err;
                        $all_done = 0;
                        last;
                    }
                    
                    if ($completed) {
                        # on instant complete, parse calls post_process itself
                        # and only returns true if that was successfull
                        $self->_complete_state($step, $state, $step_number, $pipeline, \%previous_step_outputs);
                        next;
                    }
                    else {
                        my $dispatched = $step->dispatched();
                        if (@$dispatched) {
                            foreach my $arrayref (@$dispatched) {
                                my ($cmd, $reqs, $job_args) = @$arrayref;
                                VRPipe::Submission->get(job => VRPipe::Job->get(dir => $output_root, $job_args ? (%{$job_args}) : (), cmd => $cmd), stepstate => $state, requirements => $reqs);
                            }
                        }
                        else {
                            #$self->throw("step ".$step->id." for data element ".$element->id." for pipeline setup ".$setup->id." neither completed nor dispatched anything!");
                            # it is possible for a parse to result in a different step being started over because input files were missing
                        }
                    }
                }
                
                $all_done = 0;
                last;
            }
            
            last if $incomplete_elements > $limit;
        }
        
        # capture stop
        
        return $all_done;
    }
    
    method _complete_state (VRPipe::Step $step, VRPipe::StepState $state, Int $step_number, VRPipe::Pipeline $pipeline, PreviousStepOutput $previous_step_outputs) {
        while (my ($key, $val) = each %{$step->outputs()}) {
            $previous_step_outputs->{$key}->{$step_number} = $val;
        }
        unless ($state->complete) {
            # are there a behaviours to trigger?
            my $schema = $self->result_source->schema;
            my $rs = $schema->resultset('StepBehaviour')->search({ pipeline => $pipeline->id, after_step => $step_number });
            my @behaviours;
            while (my $behaviour = $rs->next) {
                push @behaviours, $behaviour;
            }
            foreach my $behaviour (@behaviours) {
                $behaviour->behave(data_element => $state->dataelement, pipeline_setup => $state->pipelinesetup);
            }
            $state->complete(1);
            $state->update;
        }
    }
    
    #*** passing any large ref of refs to a 'method' causes a memory leak, so we
    #    must forego type checking and use plain subs for the methods that
    #    accept submission subrefs as args
    
    #method unfinished_submissions (ArrayRef[VRPipe::Submission] :$submissions?) {
    sub unfinished_submissions {
        my ($self, %args) = @_;
        my $submissions = delete $args{submissions};
        $self->throw("unexpected args") if keys %args;
        
        my @not_done;
        
        if ($submissions) {
            foreach my $sub (@$submissions) {
                push(@not_done, $sub) unless $sub->done;
            }
        }
        else {
            my $schema = $self->result_source->schema;
            my $rs = $schema->resultset('Submission')->search({ '_done' => 0 });
            while (my $sub = $rs->next) {
                push(@not_done, $sub);
            }
        }
        
        $self->debug("There are ".scalar(@not_done)." unfinished submissions to work on");
        
        return @not_done ? \@not_done : undef;
    }
    
    method handle_submissions (PositiveInt :$max_retries = 3) {
        my $submissions = $self->unfinished_submissions();
        if ($submissions) {
            $submissions = $self->check_running($submissions);
            $self->resubmit_failures(max_retries => $max_retries, submissions => $submissions);
            $self->batch_and_submit($submissions);
            return 0;
        }
        else {
            return 1;
        }
    }
    
    #method resubmit_failures (PositiveInt :$max_retries, ArrayRef[VRPipe::Submission] :$submissions) {
    sub resubmit_failures {
        my ($self, %args) = @_;
        my $max_retries = delete $args{max_retries} || $self->throw("max_retries is required");
        my $submissions = delete $args{submissions} || $self->throw("submissions is required");
        $self->throw("unexpected args") if keys %args;
        
        # this checks all the ->submissions for ones that failed, and uses the standard
        # Submission and Job methods to work out why and resubmit them as appropriate,
        # potentially with updated Requirements. max_retries defaults to 3.
        
        #*** not yet fully implemented...
        
        my $fails = 0;
        foreach my $sub (@$submissions) {
            next unless $sub->failed;
            
            if ($sub->retries >= $max_retries) {
                $self->debug("submission ".$sub->id." retried $max_retries times now, giving up");
                $fails++;
            }
            else {
                my $parser = $sub->scheduler_stdout;
                if ($parser) {
                    my $status = $parser->status;
                    #*** how can we make parser status method be present for all
                    #    scheduler outputs, and how can we make the return
                    #    values generic?
                    if ($status) {
                        if ($status eq 'MEMLIMIT') {
                            $sub->extra_memory;
                        }
                        elsif ($status eq 'RUNLIMIT') {
                            my $hrs = ceil($parser->time / 60 / 60) + 1;
                            my $current_hrs = $sub->time;
                            $sub->extra_time($hrs - $current_hrs);
                        }
                    }
                }
                $sub->retry;
            }
        }
        
        $self->debug("There were $fails permanently failed submission");
        
        return 1;
    }
    
    #method check_running (ArrayRef[VRPipe::Submission] $submissions) {
    sub check_running { 
        my ($self, $submissions) = @_;
        # this does the dance of checking if any of the currently running ->submissions
        # are approaching their time limit, and switching queues as appropriate. It
        # also checks that all running Jobs have had a recent heartbeat, and if not
        # will do the kill dance and resubmit.
        
        #*** not yet fully implemented...
        
        # update the status of each submission in case any of them finished
        my @still_not_done;
        my $c = 0;
        foreach my $sub (@$submissions) {
            my $job = $sub->job;
            $c++;
            $self->debug("loop $c, sub ".$sub->id." job ".$job->id);
            if ($job->running) {
                # user's scheduler might kill the submission if it runs too long in the
                # queue it was initially submitted to; user could do something like this
                # every 15mins:
                #if ($sub->close_to_time_limit(30)) {
                #    # where close_to_time_limit returns true if $sub->job->wall_time >
                #    # $job->requirements->time - 30.
                #    $sub->extra_time(2);
                #    # now make your scheduler recalculate the appropriate queue of a job
                #    # that takes 2 extra hours, and if the queue changes, switch the
                #    # queue:
                #    $sub->scheduler->switch_queues($sub);
                #}
                
                # check we've had a recent heartbeat
                $self->debug(" -- running...");
                if ($job->unresponsive) {
                    $self->debug(" -- unresponsive");
                    $job->kill_job;
                    $sub->update_status();
                }
            }
            elsif ($job->finished) {
                $self->debug(" -- finishing...");
                $sub->update_status();
                $self->debug(" -- success");
                next if $sub->done;
            }
            elsif ($sub->scheduled) {
                # we may be pending in the scheduler, but if we've been
                # apparently pending for ages, check the scheduler really is
                # pending us and if not reset the submission
                my $pend_time = $sub->pend_time;
                $self->debug(" -- pending for $pend_time seconds...");
                if ($pend_time > 21600) { #*** rather than check every time check_running is called after 6hrs have passed, can we only check again once every subsequent hour?
                    $self->warn("submission ".$sub->id." has been scheduled for $pend_time seconds - will try and resolve this");
                    $sub->unschedule_if_not_pending;
                }
            }
            
            push(@still_not_done, $sub);
        }
        
        $self->debug("There are ".scalar(@still_not_done)." submissions that still aren't done");
        
        return \@still_not_done;
    }
    
    #method batch_and_submit (ArrayRef[VRPipe::Submission] $submissions) {
    sub batch_and_submit {
        my ($self, $submissions) = @_;
        
        # this groups all the non-permanently-failed ->submissions into lists based
        # on shared Requirements objects, creates [PersistentArray, Requirements] for
        # each group, and uses those to do a Scheduler->submit for each group. If a
        # group consists of only 1 Submission, then it won't bother with the creating
        # a PersistentArray step.
        
        # we can have step-specific limits on how many subs to run at once; see
        # if any are in place and set them globally for all setups
        my %step_limits;
        my $do_step_limits = 0;
        foreach my $setup ($self->setups) {
            if ($setup->active) {
                my $user_opts = $setup->options;
                
                my @stepms = $setup->pipeline->step_members;
                foreach my $stepm (@stepms) {
                    my $step = $stepm->step;
                    my $step_name = $step->name;
                    
                    my $limit = $user_opts->{$step_name.'_max_simultaneous'};
                    unless ($limit) {
                        $limit = $step->max_simultaneous;
                    }
                    
                    if ($limit) {
                        if (! defined $step_limits{$step_name} || $limit < $step_limits{$step_name}) {
                            $step_limits{$step_name} = $limit;
                        }
                        $do_step_limits = 1;
                    }
                }
            }
        }
        
        # if any step limits are in place we first have to see how many of each
        # step are currently in the scheduler
        my %step_counts;
        if ($do_step_limits) {
            my $schema = $self->result_source->schema;
            my $rs = $schema->resultset('Submission')->search({ '_done' => 0, '_sid' => { '!=', undef } });
            while (my $sub = $rs->next) {
                $step_counts{$sub->stepstate->stepmember->step->name}++;
            }
        }
        
        my %batches;
        my $limit = $self->global_limit;
        my $count = 0;
        foreach my $sub (@$submissions) {
            next if ($sub->sid || $sub->done || $sub->failed);
            
            # step limit handling
            if ($do_step_limits) {
                my $step_name = $sub->stepstate->stepmember->step->name;
                if (exists $step_limits{$step_name}) {
                    $step_counts{$step_name}++;
                    next if $step_counts{$step_name} > $step_limits{$step_name};
                }
            }
            
            # global limit handling
            $count++;
            last if $count >= $limit;
            
            # we're good to batch this one
            push(@{$batches{$sub->requirements->id}}, $sub);
        }
        
        my $scheduler = VRPipe::Scheduler->get;
        while (my ($req_id, $subs) = each %batches) {
            my $sid = $scheduler->submit(array => $subs,
                                         requirements => VRPipe::Requirements->get(id => $req_id));
        }
    }
}

1;