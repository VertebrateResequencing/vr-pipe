
=head1 NAME

VRPipe::Manager - methods for managing the execution of pipelines

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This is the main module used to discover what work needs to be done ('trigger'
each L<VRPipe::PipelineSetup>) and then to 'dispatch' that work out to the
system's job scheduler to actually run the command lines on the compute
cluster.

The B<vrpipe-trigger_pipelines> and B<vrpipe-dispatch_pipelines> scripts call
methods of this module.

Note that the current implementation is slow and inefficient, with lots of
serial looping over thousands of objects. A radical overhaul is planned.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::Manager extends VRPipe::Persistent {
    use Parallel::ForkManager;
    use Sys::CPU;
    use POSIX qw(ceil);
    
    our $DEFAULT_MAX_PROCESSES = Sys::CPU::cpu_count();
    our %step_limits;
    our %setups_with_step_limits;
    our %setups_checked_for_step_limits;
    our %step_limit_instigators;
    our $do_step_limits = 0;
    
    has '_running' => (is      => 'rw',
                       isa     => 'Bool',
                       traits  => ['VRPipe::Persistent::Attributes'],
                       default => 0);
    
    has '_run_start' => (is          => 'rw',
                         isa         => Datetime,
                         traits      => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    has '_instance_run_start' => (is  => 'rw',
                                  isa => Datetime);
    
    has 'global_limit' => (is  => 'rw',
                           isa => PositiveInt);
    
    has '_created_pipelines' => (is  => 'rw',
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
        return $self->create(global_limit => $global_limit);
    }
    
    around create (ClassName|Object $self: Persistent :$id?, PositiveInt :$global_limit = 500) {
        # our db-storage needs are class-wide, so we only have 1 row in our
        # table
        return $self->$orig(id => 1, global_limit => $global_limit);
    }
    
    method setups (Str :$pipeline_name?) {
        my @setups;
        foreach my $ps (VRPipe::PipelineSetup->search({}, { prefetch => 'pipeline' })) {
            my $p = $ps->pipeline;
            if ($pipeline_name) {
                next unless $p->name eq $pipeline_name;
            }
            
            # before spooling, whilst we are still in a single process, make
            # sure that all pipelines have had their step members created using
            # the steps method. *** steps() currently gives a memory leak, so
            # that's the other reason we are careful to only call it once per
            # pipeline
            my $ps_id   = $ps->id;
            my $created = $self->_created_pipelines;
            unless (exists $created->{ $p->id }) {
                $p->steps;
                $created->{ $p->id } = 1;
                $self->_created_pipelines($created);
            }
            
            # we can have step-specific limits on how many subs to run at
            # once; see if any are in place and set them globally for all
            # setups
            if ($ps->active && !exists $setups_checked_for_step_limits{$ps_id}) {
                my $user_opts = $ps->options;
                
                my @stepms = $ps->pipeline->step_members;
                foreach my $stepm (@stepms) {
                    my $step      = $stepm->step;
                    my $step_name = $step->name;
                    
                    my $limit = $user_opts->{ $step_name . '_max_simultaneous' };
                    unless ($limit) {
                        $limit = $step->max_simultaneous;
                    }
                    
                    if ($limit) {
                        if (!defined $step_limits{$step_name} || $limit <= $step_limits{$step_name}) {
                            $step_limits{$step_name}                                = $limit;
                            $setups_with_step_limits{$ps_id}->{$step_name}          = $limit;
                            $step_limit_instigators{$step_name}->{$limit}->{$ps_id} = 1;
                            $do_step_limits                                         = 1;
                        }
                    }
                }
                
                $setups_checked_for_step_limits{$ps_id} = 1;
            }
            elsif (!$ps->active && exists $setups_with_step_limits{$ps_id}) {
                my @limited_steps = keys %{ $setups_with_step_limits{$ps_id} };
                while (my ($step_name, $limit) = each %{ $setups_with_step_limits{$ps_id} }) {
                    delete $step_limit_instigators{$step_name}->{$limit}->{$ps_id};
                    my $current_limit = $step_limits{$step_name};
                    unless (keys %{ $step_limit_instigators{$step_name}->{$current_limit} }) {
                        delete $step_limit_instigators{$step_name}->{$current_limit};
                        my @other_limits = sort { $a <=> $b } keys %{ $step_limit_instigators{$step_name} };
                        if (@other_limits) {
                            $step_limits{$step_name} = $other_limits[0];
                        }
                        else {
                            delete $step_limits{$step_name};
                            delete $step_limit_instigators{$step_name};
                        }
                    }
                }
                
                delete $setups_with_step_limits{$ps_id};
                delete $setups_checked_for_step_limits{$ps_id};
                
                unless (keys %setups_with_step_limits) {
                    $do_step_limits = 0;
                }
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
        $fm->run_on_finish(
            sub {
                my $data_structure_reference = $_[5];
                if (defined($data_structure_reference)) {
                    $retrieved_responses{$$data_structure_reference} = 1;
                }
                else {
                    $self->throw("spool child failed to return anything");
                }
            });
        
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
        
        my $setup_id     = $setup->id;
        my $pipeline     = $setup->pipeline;
        my @step_members = $pipeline->step_members;
        my $num_steps    = scalar(@step_members);
        
        my $datasource  = $setup->datasource;
        my $output_root = $setup->output_root;
        $self->make_path($output_root);
        my $all_done            = 1;
        my $incomplete_elements = 0;
        my $limit               = $self->global_limit;
        my %step_counts;
        
        # We don't pass a limit to incomplete_element_states since that could
        # return all estates where all their submissions have failed 3 times.
        # To avoid doing wasted cycles of the estates loop we do make use of
        # limit though. The real limiting happens in handle_submissions()
        my $pager = $datasource->incomplete_element_states($setup);
        
        $self->debug("Spool for setup " . $setup->id . " got " . $pager->total_entries . " incomplete element states");
      
      ELOOP: while (my $estates = $pager->next) {
            foreach my $estate (@$estates) {
                my $element         = $estate->dataelement;
                my $completed_steps = $estate->completed_steps;
                $self->debug("setup $setup_id | incomplete element loop $incomplete_elements for element id " . $element->id . " which has completed $completed_steps steps of $num_steps");
                next if $completed_steps == $num_steps; # (rare, if ever?)
                
                $incomplete_elements++;
                
                my %previous_step_outputs;
                my $already_completed_steps = 0;
                foreach my $member (@step_members) {
                    my $step_number = $member->step_number;
                    my $state = VRPipe::StepState->create(stepmember    => $member,
                                                          dataelement   => $element,
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
                    
                    my $step_name      = $step->name;
                    my $inc_step_count = 0;
                    if (exists $step_limits{$step_name}) {
                        $inc_step_count = 1;
                        if ($step_counts{$step_name} && $step_counts{$step_name} > $step_limits{$step_name}) {
                            $all_done = 0;
                            last;
                        }
                    }
                    
                    # have we previously done the dispatch dance and are currently
                    # waiting on submissions to complete?
                    my @submissions = $state->submissions;
                    if (@submissions) {
                        my $unfinished = VRPipe::Submission->search({ '_done' => 0, stepstate => $state->id });
                        unless ($unfinished) {
                            my $ok = $step->post_process();
                            if ($ok) {
                                $self->debug(" we just completed all the submissions from a previous parse");
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
                            $self->debug(" we have $unfinished unfinished submissions from a previous parse");
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
                            else {
                                $step_counts{$step_name}++ if $inc_step_count;
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
                            $self->debug(" parsing the step resulted in instant completion");
                            next;
                        }
                        else {
                            my $dispatched = $step->dispatched();
                            if (@$dispatched) {
                                foreach my $arrayref (@$dispatched) {
                                    my ($cmd, $reqs, $job_args) = @$arrayref;
                                    my $sub = VRPipe::Submission->create(job => VRPipe::Job->create(dir => $output_root, $job_args ? (%{$job_args}) : (), cmd => $cmd), stepstate => $state, requirements => $reqs);
                                    $self->debug(" parsing the step made new submission " . $sub->id . " with job " . $sub->job->id);
                                }
                                $step_counts{$step_name}++ if $inc_step_count;
                            }
                            else {
                                $self->debug("step " . $step->id . " for data element " . $element->id . " for pipeline setup " . $setup->id . " neither completed nor dispatched anything!");
                                # it is possible for a parse to result in a different step being started over because input files were missing
                            }
                        }
                    }
                    
                    $all_done = 0;
                    last;
                }
                
                last ELOOP if $incomplete_elements > $limit;
            }
        }
        
        # capture stop
        
        return $all_done;
    }
    
    method _complete_state (VRPipe::Step $step, VRPipe::StepState $state, Int $step_number, VRPipe::Pipeline $pipeline, PreviousStepOutput $previous_step_outputs) {
        while (my ($key, $val) = each %{ $step->outputs() }) {
            $previous_step_outputs->{$key}->{$step_number} = $val;
        }
        unless ($state->complete) {
            # are there a behaviours to trigger?
            foreach my $behaviour (VRPipe::StepBehaviour->search({ pipeline => $pipeline->id, after_step => $step_number })) {
                $behaviour->behave(data_element => $state->dataelement, pipeline_setup => $state->pipelinesetup);
            }
            
            # add to the StepStats
            foreach my $submission ($state->submissions) {
                my $sched_stdout = $submission->scheduler_stdout || next;
                my $memory = ceil($sched_stdout->memory || $submission->memory);
                my $time   = ceil($sched_stdout->time   || $submission->time);
                VRPipe::StepStats->create(step => $step, pipelinesetup => $state->pipelinesetup, submission => $submission, memory => $memory, time => $time);
            }
            
            $state->complete(1);
            $state->update;
        }
    }
    
    method handle_submissions (Int :$max_retries = 3) {
        # call setups to find step limits
        $self->setups;
        
        my $unfinished = VRPipe::Submission->search({ '_done' => 0 }) || 0;
        $self->debug("There are $unfinished unfinished submissions to work on");
        
        if ($unfinished) {
            $self->check_running();
            $self->resubmit_failures(max_retries => $max_retries);
            $self->batch_and_submit();
            return 0;
        }
        else {
            return 1;
        }
    }
    
    method check_running {
        # update the status of each submission in case any of them finished
        my $pager = VRPipe::Submission->search_paged({ '_done' => 0, '_failed' => 0 }, { prefetch => 'job' });
        
        my $still_not_done = 0;
        my $c              = 0;
        while (my $subs = $pager->next) {
            foreach my $sub (@$subs) {
                my $job = $sub->job;
                $c++;
                $self->debug("loop $c, sub " . $sub->id . " job " . $job->id);
                if ($job->running) {
                    # check we've had a recent heartbeat
                    $self->debug(" -- running...");
                    if ($job->unresponsive) {
                        $self->debug(" -- unresponsive");
                        $job->kill_job;
                        $sub->update_status();
                    }
                    elsif ($sub->close_to_time_limit(30)) {
                        # user's scheduler might kill the submission if it runs too
                        # long in the queue it was initially submitted to; avoid
                        # this by changing queue as necessary
                        $self->debug(" -- within 30mins of time limit, will increase by 2hrs...");
                        $sub->extra_time(2);
                    }
                }
                elsif ($job->finished) {
                    $self->debug(" -- finishing by calling sub->update_status...");
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
                        $self->warn("submission " . $sub->id . " has been scheduled for $pend_time seconds - will try and resolve this");
                        $sub->unschedule_if_not_pending;
                    }
                }
                
                $still_not_done++;
            }
        }
        
        $self->debug("There are $still_not_done submissions that still aren't done");
    }
    
    method resubmit_failures (Int :$max_retries) {
        # this checks all failed submissions, and uses the standard Submission
        # and Job methods to work out why and resubmit them as appropriate,
        # potentially with updated Requirements. max_retries defaults to 3.
        
        my $pager = VRPipe::Submission->search_paged({ _failed => 1, retries => { '<' => $max_retries } }, { prefetch => [qw(scheduler requirements)] });
        
        while (my $subs = $pager->next) {
            foreach my $sub (@$subs) {
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
                            my $hrs         = ceil($parser->time / 60 / 60) + 1;
                            my $current_hrs = $sub->time;
                            $sub->extra_time($hrs - $current_hrs);
                        }
                    }
                }
                $sub->retry;
            }
        }
    }
    
    method batch_and_submit {
        # this groups all the non-permanently-failed ->submissions into lists
        # based on shared Requirements objects, creates [PersistentArray,
        # Requirements] for each group, and uses those to do a Scheduler->submit
        # for each group. If a group consists of only 1 Submission, then it
        # won't bother with the creating a PersistentArray step.
        
        # if any step limits are in place we first have to see how many of each
        # step are currently in the scheduler
        my %step_counts;
        if ($do_step_limits) {
            my $array_ref = VRPipe::Submission->get_column_values(['step.name', 'count(step.name)'],
                                                                  { '_done' => 0, '_failed' => 0, '_sid' => { '!=', undef } },
                                                                  {  join     => { stepstate => { stepmember => 'step' } },
                                                                     group_by => ['step.name'] });
            foreach my $vals (@$array_ref) {
                $step_counts{ $vals->[0] } = $vals->[1];
            }
        }
        
        # we want to limit the number of jobs we have running based on the
        # global limit, so we'll count the number of submissions we batch and
        # end the submission loop when over the limit. We'll also start the
        # count at the number of currently running jobs.
        my $limit = $self->global_limit;
        my $count = VRPipe::Job->search({ 'running' => 1 });
        
        # we want to now find all submissions that are not done or failed and
        # that haven't already been submitted. Because another process could
        # create new submissions which would change our results during the loop,
        # we avoid unecessary (possibly infinite) loop restarts by Pager by
        # first getting the most recently created submission and then searching
        # for ids less than that.
        my ($last_sub_id) = VRPipe::Submission->get_column_values('id', {}, { order_by => { -desc => 'id' }, rows => 1 });
        my $pager = VRPipe::Submission->search_paged({ '_done' => 0, '_failed' => 0, '_sid' => undef, 'me.id' => { '<=' => $last_sub_id } }, { order_by => 'requirements', prefetch => [qw(job requirements)] }, 1000);
        
        my $scheduler = VRPipe::Scheduler->get;
      SLOOP: while (my $subs = $pager->next) {
            my %batches;
            foreach my $sub (@$subs) {
                # if the job is a block_and_skip_if_ok job, we don't actually
                # block because of race condition issues, and because of fail
                # and restart issues. Instead we always only actually submit if
                # this submission is the first submission created for this $job.
                # This also saves us creating lots of submissions for a single
                # job, which is confusing and wasteful.
                my $job = $sub->job;
                if ($job->block_and_skip_if_ok) {
                    my ($first_sub) = VRPipe::Submission->search({ 'job' => $job->id }, { rows => 1, order_by => { -asc => 'id' } });
                    next unless $first_sub->id == $sub->id;
                }
                
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
                last SLOOP if $count >= $limit;
                
                # we're good to batch this one
                push(@{ $batches{ $sub->requirements->id } }, $sub);
            }
            
            while (my ($req_id, $subs) = each %batches) {
                $self->debug("_~_ Batched and submitted an array of " . scalar(@$subs) . " subs");
                my $sid = $scheduler->submit(array        => $subs,
                                             requirements => $subs->[0]->requirements);
            }
        }
    }
}

1;
