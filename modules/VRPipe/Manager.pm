use VRPipe::Base;

class VRPipe::Manager extends VRPipe::Persistent {
    use Parallel::ForkManager;
    use Sys::CPU;
    
    our $DEFAULT_MAX_PROCESSES = Sys::CPU::cpu_count();
    
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
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
    
    # public getters for our private attributes
    method running {
        return $self->_running;
    }
    method run_start {
        return $self->_run_start;
    }
    
    __PACKAGE__->make_persistent();
    
    around get (ClassName|Object $self: Persistent :$id?) {
        # our db-storage needs are class-wide, so we only have 1 row in our
        # table
        return $self->$orig(id => 1);
    }
    
    method setups (Str :$pipeline_name?) {
        my $rs = $self->result_source->schema->resultset('PipelineSetup');
        my @setups;
        while (my $ps = $rs->next) {
            if ($pipeline_name) {
                my $p = $ps->pipeline;
                next unless $p->name eq $pipeline_name;
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
        
        my $pipeline = $setup->pipeline;
        my @step_members = $pipeline->steps;
        
        my $datasource = $setup->datasource;
        my $output_root = $setup->output_root;
        my $all_done = 1;
        foreach my $element ($datasource->elements) {
            my %previous_step_outputs;
            foreach my $member (@step_members) {
                my $state = VRPipe::StepState->get(stepmember => $member,
                                                   dataelement => $element,
                                                   pipelinesetup => $setup);
                
                my $step = $member->step(previous_step_outputs => \%previous_step_outputs, step_state => $state);
                if ($state->complete) {
                    $self->_complete_state($step, $state, \%previous_step_outputs);
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
                            $self->_complete_state($step, $state, \%previous_step_outputs);
                            next;
                        }
                        else {
                            $self->throw("submissions completed, but post_process failed");
                        }
                    }
                }
                else {
                    # this is the first time we're looking at this step for
                    # this data member for this pipelinesetup
                    my $completed = $step->parse();
                    if ($completed) {
                        # on instant complete, parse calls post_process itself
                        # and only returns true if that was successfull
                        $self->_complete_state($step, $state, \%previous_step_outputs);
                        next;
                    }
                    else {
                        my $dispatched = $step->dispatched();
                        if (@$dispatched) {
                            foreach my $arrayref (@$dispatched) {
                                my ($cmd, $reqs) = @$arrayref;
                                my $sub = VRPipe::Submission->get(job => VRPipe::Job->get(cmd => $cmd, dir => $output_root), stepstate => $state, requirements => $reqs);
                            }
                        }
                        else {
                            $self->throw("step ".$step->id." for data element ".$element->id." for pipeline setup ".$setup->id." neither completed nor dispatched anything!");
                        }
                    }
                }
                
                $all_done = 0;
                last;
            }
        }
        
        # capture stop
        
        return $all_done;
    }
    
    method _complete_state (VRPipe::Step $step, VRPipe::StepState $state, HashRef $previous_step_outputs) {
        while (my ($key, $val) = each %{$step->outputs()}) {
            $previous_step_outputs->{$key} = $val;
        }
        unless ($state->complete) {
            $state->complete(1);
            $state->update;
        }
    }
    
    method unfinished_submissions (ArrayRef[VRPipe::Submission] :$submissions?) {
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
    
    method resubmit_failures (PositiveInt :$max_retries, ArrayRef[VRPipe::Submission] :$submissions) {
        # this checks all the ->submissions for ones that failed, and uses the standard
        # Submission and Job methods to work out why and resubmit them as appropriate,
        # potentially with updated Requirements. max_retries defaults to 3.
        
        #*** not yet implemented...
        
        return 1;
    }
    
    method check_running (ArrayRef[VRPipe::Submission] $submissions) {
        # this does the dance of checking if any of the currently running ->submissions
        # are approaching their time limit, and switching queues as appropriate. It
        # also checks that all running Jobs have had a recent heartbeat, and if not
        # will do the kill dance and resubmit.
        
        #*** not yet fully implemented...
        
        # update the status of each submission in case any of them finished
        my @still_not_done;
        foreach my $sub (@$submissions) {
            if ($sub->job->running) {
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
            }
            elsif ($sub->job->finished) {
                $sub->update_status();
                next;
            }
            else {
                # we must be pending in the scheduler
                #my $pend_time = $sub->pend_time;
                # this also works when we're running/finished by taking the difference
                # between $job->start_time and $self->scheduled. Otherwise its the
                # difference between the time now and $self->scheduled.
            }
            
            push(@still_not_done, $sub);
        }
        
        return \@still_not_done;
    }
    
    method batch_and_submit (ArrayRef[VRPipe::Submission] $submissions) {
        # this groups all the non-permanently-failed ->submissions into lists based
        # on shared Requirements objects, creates [PersistentArray, Requirements] for
        # each group, and uses those to do a Scheduler->submit for each group. If a
        # group consists of only 1 Submission, then it won't bother with the creating
        # a PersistentArray step.
        
        my %batches;
        foreach my $sub (@$submissions) {
            next if ($sub->sid || $sub->done || $sub->failed);
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