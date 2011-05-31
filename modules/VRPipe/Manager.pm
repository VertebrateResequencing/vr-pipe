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
        foreach my $setup (@$setups) {
            if ($setup->active) {
                $fm->start and next; # fork
                
                my $num_of_jobs_still_to_go = $self->spool($setup);
                
                $fm->finish; # exit in the child process
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
        
        warn "spooling setup ", $setup->id, "\n";
        my $pipeline = $setup->pipeline;
        my @step_members = $pipeline->steps;
        warn "got steps (@step_members)\n";
        
        my $datasource = $setup->datasource;
        my $output_root = $setup->output_root;
        foreach my $element ($datasource->elements) {
            warn "got element $element with id ", $element->id, "\n";
            my $input = $element;
            foreach my $member (@step_members) {
                my $state = VRPipe::StepState->get(stepmember => $member,
                                                   dataelement => $element,
                                                   pipelinesetup => $setup);
                
                warn "got state $state with id ", $state->id, "\n";
                
                my $step = $member->step(input => $input, step_state => $state);
                if ($state->complete) {
                    $input = $step->outputs();
                    next;
                }
                
                # have we previously done the dispatch dance and are currently
                # waiting on submissions to complete?
                if (0) {
                    # *** ???
                }
                else {
                    # this is the first time we're looking at this step for
                    # this data member for this pipelinesetup
                    my $completed = $step->parse();
                    if ($completed) {
                        # on instant complete, parse calls post_process itself
                        # and only returns true if that was successfull
                        $state->complete(1);
                        $state->update;
                        next;
                    }
                    else {
                        my @dispatched = $step->dispatched();
                        if (@dispatched) {
                            # this returns a list of [$cmd_line_string, $requirements_object] refs that
                            # you could use to create VRPipe::JobManager::Submission objects and
                            # eventually actually get the cmd_line to run.
                            
                            
                            # You'll keep records on how those Submissions do, and when
                            # they're all successfull you'll do:
                            #$step->post_process();
                        }
                        else {
                            $self->throw("step ".$step->id." for data element ".$element->id." for pipeline setup ".$setup->id." neither completed nor dispatched anything!");
                        }
                    }
                }
                
                last;
            }
        }
        
        # capture stop
    }
}

1;

=pod

 # get a list of all VRPipe::PipelineManager::Setup objects:
 my @setups = $man->setups;
 # limited to those for a certain pipeline module:
 @setups = $man->setups(module => 'VRPipe::Pipelines::Mapping');

 # trigger/find out about the pipelines:
 foreach my $set (@setups) {
    if ($set->active) {
        # run the pipeline (if it is complete, nothing happens; since new data
        # could have turned up, there is no general concept of a pipeline ever
        # automatically being considered 'complete' so there is no "first check
        # we haven't already completed" method at this level):
        my $num_of_jobs_still_to_go = $man->run($set);
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
        #
        # If run($pip, stdout => 1) is called, then it outputs this kind of
        # stuff in a pretty way to STDOUT as it goes along. Otherwise this info
        # is queryable after run($pip) returns:
        my %status = $man->status($set);
        
        # calling run() should be very cheap, since for the most part it will
        # just be doing lookups into the Persistent db where our state is
        # stored. However, if a step method has to actually run, then it may
        # not be fast enough for constructing a live-updated status webpage. For
        # that purpose, ->status($set) can be called without calling ->run($set)
        # first, in which case it will get all the stats it can by going through
        # the same loops as run(), except that it won't call ->run_step (so
        # will be able to say that a step has not yet started, but not how
        # many jobs there are to go for that step). ->status($set,
        # overall_jobs_only => 1) takes a shortcut and just directly grabs the
        # Submission objects from $set->children and generates a simple set of
        # completed/running/to-be-run/failed overall stats. NB: we never know
        # about steps that we haven't gotten to, so the total number of jobs
        # is never accurate until the last step has been run on all data
        # elements.
        
        # a simple hash from status isn't powerful/easy enough to allow things
        # like showing failed jobs and error messages in a web-frontend. For
        # that purpose, we have a set of methods that would let a user drill
        # down and discover what's going wrong with their pipeline:
        my @failed_steps = $man->failed_steps($set);
        # a 'failed' step is one that has submissions that have finished but
        # not ok.
        foreach my $step (@failed_steps) {
            # meta-information (like the dataelement_key) can be extracted in
            # the normal way from a step for display
            
            my %outputs = $man->failed_outputs($step);
            # this gets the last STDOUT and STDERR from every failed Submission
            # associated with $step and returns a hash with Submission->id
            # keys and {stdout => 'string', stderr => 'string'} values.
            
            # you could imagine that a web-frontend would have a button that
            # called:
            $step->reset;
            # that the user might use if they looked at the outputs and fixed
            # the problem.
        }
        
        # perhaps the user fixed a problem and now wants to reset all the
        # steps that failed in one go:
        $man->reset(@failed_steps);
        
        # perhaps something really stupid and wrong happened with a pipeline and
        # you just want to start everything over completely from scratch:
        $man->reset($set);
        # this grabs all the Step objects from $set->children and does ->reset
        # on them.
        
        # run-time/memory-usage summary stats can be found on a per-step
        # basis, averaged over every pipeline that ran that step, and also
        # per pipeline:
        my %stats = $man->step_stats(); # averaged over all pipelines
        %stats = $man->step_stats($set); # for just this pipeline config
        # %stats == (step_name => { walltime => $seconds, memory => $mb });
    }
    elsif (time() - $set->last_updated > 7776000) {
        # it's been inactive for ~3months; perhaps we want to trash old stuff?
        $man->expunge($set);
        # expunges all Persistent entries in the Peristent db that were created
        # during the course of this pipeline setup, that are not used by another
        # pipeline setup
    }
 }

 # Having ->run at least 1 pipeline setup, we'll have a bunch of Submission
 # objects in the Persistent db, but none of them will have been submitted. We
 # deal with submitting, coping with failed Submissions, retrying Submissions
 # etc. in the following few methods, which don't need to know anything about
 # the pipelines that spawned them. This is completly generic, and you could
 # imagine running a script that calls these methods every 10mins in a cron job
 # or something.

 # get all Submissions that have not completed ok (used internally by the
 # subsequent methods, answer cached once per instance of Manager):
 my @incomplete_submissions = $man->submissions();
 # again, this is all incomplete Submission objects, regardless of what pipeline
 # or step or process spawned them. submissions() will block if the running()
 # boolean is true, to give a run() time to finish, so we'll see all of its
 # Submission objects. After a reasonable amount of time it will stop blocking
 # though, incase a run() went bad and failed to turn off running().

 # deal with failed Submissions:
 $man->resubmit_failures(max_retries => 3);
 # this checks all the ->submissions for ones that failed, and uses the standard
 # Submission and Job methods to work out why and resubmit them as appropriate,
 # potentially with updated Requirements. max_retries defaults to 3.

 $man->check_running;
 # this does the dance of checking if any of the currently running ->submissions
 # are approaching their time limit, and switching queues as appropriate. It
 # also checks that all running Jobs have had a recent heartbeat, and if not
 # will do the kill dance and resubmit.

 # batch Submissions into arrays and submit them:
 $man->submit;
 # this groups all the non-permanently-failed ->submissions into lists based
 # on shared Requirements objects, creates [PersistentArray, Requirements] for
 # each group, and uses those to do a Scheduler->submit for each group. If a
 # group consists of only 1 Submission, then it won't bother with the creating
 # a PersistentArray step. Suppling ->submit($percent) makes it consider only
 # that percentage of incomplete Submissions. This would be useful in a script
 # that would be run by, say, 4 team members: the script would say ->submit(25),
 # and the workload would be shared amongst the members, helpful when dealing
 # with fair-share systems in schedulers.
 
 # expunge all the PersistentArray objects that submit() has made in the past,
 # but are now defunct (they reference Submissions that have all finished OK):
 $man->cleanup;

=cut