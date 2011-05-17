use VRPipe::Base;

class VRPipe::Scheduler extends VRPipe::Persistent {
    use VRPipe::Config;
    my $vrp_config = VRPipe::Config->new();
    use VRPipe::Persistent::SchemaBase;
    
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'type' => (is => 'rw',
                   isa => Varchar[64],
                   builder => 'default_type',
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1,
                   allow_key_to_default => 1);
    
    has 'output_root' => (is => 'rw',
                          isa => Varchar[64],
                          builder => 'default_output_root',
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1,
                          allow_key_to_default => 1);
    
    method default_type (ClassName|Object $self:) {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment.'_scheduler';
        return $vrp_config->$method_name();
    }
    
    method default_output_root (ClassName|Object $self:) {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment.'_scheduler_output_root';
        return $vrp_config->$method_name();
    }
    
=head2 submit

 Title   : submit
 Usage   : my $scheduled_id = $obj->submit(submission => $VRPipe_submission);
           $scheduled_id = $obj->submit(array => $VRPipe_persistentarray,
                                        requirements => $VRPipe_requirements);
 Function: Submit one or more jobs to the scheduler (queue them to be run).
 Returns : int (id from the scheduler)
 Args    : submission => VRPipe::Submission instance
           -OR-
           array => VRPipe::PersistentArray | [VRPipe::Submission instances]
           requirements => VRPipe::Requirements instance

=cut
    method submit (VRPipe::Submission :$submission?, VRPipe::Requirements :$requirements?, PersistentArray|ArrayRefOfPersistent :$array?) {
        #my $scheduled_id = $sch->submit(submission => $VRPipe_jobmanager_submission);
        # this gets requirements from the Submission object. It also auto-sets
        # $submission->sid() to $scheduled_id (the return value, which is the value
        # returned by the scheduler). _hid gets set to submission id and _aid
        # to 0.
        
        #$scheduled_id = $sch->submit(array => $VRPipe_persistentarray,
        #                      requirements => $VRPipe_jobmanager_requirements);
        # for running more than one job in an array, you pass a PersistentArray object
        # and a Requirments object (ie. where you have arranged that all the Submission
        # objects in the array share this same Requirements). It gets each Submission
        # object from the PersistentArray and auto-sets $submission->sid() to the
        # $scheduled_id with the PersistentArray id and index in _hid and _aid.
        # In both cases, it claims all Submission objects prior to interacting with
        # the scheduler, and if the scheduler has a problem and the submission fails,
        # we release all the Submission objects.
        
        my @submissions;
        my $for;
        my $aid = -1;
        if ($submission) {
            push(@submissions, $submission);
            $requirements ||= $submission->requirements;
            $for = $submission;
        }
        elsif ($array) {
            if (ref($array) eq 'ARRAY') {
                $array = VRPipe::PersistentArray->get(members => $array);
            }
            push(@submissions, @{$array->members});
            $for = $array;
            $aid++;
        }
        $self->throw("at least one submission and requirements must be supplied") unless @submissions && $requirements;
        
        # generate a command line that will submit to the scheduler
        my $cmd_line = $self->build_command_line(requirements => $requirements, for => $for);
        
        # claim all submission objects, associating them with the hashing id,
        # then attempt the submit and set their sid on success or release on
        # failure
        my $schema = $self->result_source->schema;
        my $sid;
        try {
            $sid = $schema->txn_do(sub {
                my $all_claimed = 1;
                foreach my $sub (@submissions) {
                    my $claimed = $sub->claim;
                    unless ($claimed) {
                        $all_claimed = 0;
                        last;
                    }
                    $sub->_hid($for->id);
                    $sub->_aid(++$aid);
                }
                
                if ($all_claimed) {
                    foreach my $sub (@submissions) {
                        $sub->update;
                    }
                    
                    my $got_sid = $self->get_sid($cmd_line);
                    
                    if ($got_sid) {
                        foreach my $sub (@submissions) {
                            $sub->sid($got_sid);
                        }
                        return $got_sid;
                    }
                    else {
                        foreach my $sub (@submissions) {
                            $sub->release;
                        }
                        die "failed to submit to scheduler";
                    }
                }
                else {
                    die "failed to claim all submissions";
                }
            });
        }
        catch ($err) {
            $self->throw("Rollback failed!") if ($err =~ /Rollback failed/);
            $self->throw("Failed to claim & submit: $err");
        }
        
        return $sid;
        
        # When dealing with a PersistentArray, each Submission will get its own output
        # files in the output_dir() based on which index it was in the PersistentArray.
        # You access a particular Submission's output files with:
        #my $o_file = $sch->stdout($VRPipe_jobmanager_submission);
        # and likewise with ->stderr. It works out the hashing-id and index from the
        # Submission->sid (for PersistentArray submits) or Submission->id (for single
        # Submission submits). Having dealt with the output files, you can delete
        # them:
        #$sch->unlink_output($VRPipe_jobmanager_submission);
        # this also prunes empy dirs.
        
        # the command that submit() actually schedules in the scheduler is a little
        # perl -e that takes the referenced Submission id (working it out from the
        # PersistentArray object and index if this was an PersistentArray submit),
        # pulls out the Job and does ->run on that.
        
        # though some schedulers may have some kind of limit in place on the maximum
        # number of jobs a user or the system as a whole can keep scheduled at once,
        # we deliberatly don't model or support that limit here. If we do go over the
        # limit we treat it as any other error from the scheduler: ->submit would just
        # return false and ->last_submit_error would give the error as a string. The
        # user is then free to try the submit again later. To help the user decide if
        # we're under the limit or not before attempting a submit:
        #my $max_jobs = $sch->max_jobs;
        # this is the get-only max jobs as set in SiteConfig
        #my $current_jobs = $sch->queued_jobs;
        # this is the total number of jobs being tracked in the scheduling system for
        # the current user.
        # For testing purposes, VRPipe::Schedulers::Local will accept
        # a practically-infinite-sized list (max_jobs is v.high), and handle tracking
        # and serially running the jobs one-at-a-time itself.
    }
    
    method get_sid (Str $cmd) {
        #*** supposed to be implemented in eg. VRPipe::Schedulers::LSF if
        #    $self->type eq 'LSF'; hard-coded to something LSF&Sanger specific
        #    for now
        my $output = `$cmd`;
        my ($sid) = $output =~ /Job \<(\d+)\> is submitted/;
        if ($sid) {
            return $sid;
        }
        else {
            $self->throw("Failed to submit to scheduler");
        }
    }
    
    method build_command_line (VRPipe::Requirements :$requirements, PersistentObject :$for) {
        # figure out where STDOUT & STDERR of the scheduler should go
        my $output_dir = $self->output_dir($for);
        
        # the command we want each node to execute once it is run by the
        # scheduler
        my $self_id = $self->id;
        my $node_run_args = $for->isa('VRPipe::PersistentArray') ? "index => shift, array => " : "submission => ";
        $node_run_args .= $for->id;
        my $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
        my $cmd = qq[perl -MVRPipe::Persistent::SchemaBase -MVRPipe::Scheduler -e "VRPipe::Persistent::SchemaBase->database_deployment(q[$deployment]); VRPipe::Scheduler->get(id => $self_id)->run_on_node($node_run_args);"];
        
        return join(' ', $self->submit_command, $self->submit_args(requirements => $requirements, output_dir => $output_dir, cmd => $cmd, multiple => $for->isa('VRPipe::PersistentArray')));
    }
    
    method output_dir (PersistentObject $for) {
        my $root_dir = $self->output_root;
        my $hashing_string = ref($for).'::'.$for->id;
        
        #*** supposed to be placed in subdirs based on a hashing of
        #    $hashing_string, but for now just stick it in a subdir directly
        #    named after $hashing_string
        $hashing_string =~ s/\:/_/g;
        my $output_dir = dir($root_dir, $hashing_string);
        
        $output_dir->mkpath;
        return $output_dir;
    }
    
    method submit_command {
        #*** supposed to be implemented in eg. VRPipe::Schedulers::LSF if
        #    $self->type eq 'LSF'; hard-coded to something LSF&Sanger specific
        #    for now
        return 'bsub';
    }
    
    method submit_args (VRPipe::Requirements :$requirements, Dir :$output_dir, Str :$cmd, Bool :$multiple) {
        #*** supposed to be implemented in eg. VRPipe::Schedulers::LSF if
        #    $self->type eq 'LSF'; hard-coded to something LSF&Sanger specific
        #    for now
        
        # access the requirments object and build up the string based on memory,
        # time, cpu etc.
        my $queue = $self->determine_queue($requirements);
        # *** ...
        my $requirments_string = "-q $queue";
        
        # work out the scheduler output locations and how to pass on the
        # scheduler array index to the perl cmd
        my $index_spec = '';
        my $ofile = file($output_dir, 'stdout');
        my $efile = file($output_dir, 'stderr');
        my $output_string;
        if ($multiple) {
            $index_spec = ''; #*** something that gives the index to be shifted into perl -e
            $output_string = "-o $ofile -e $efile"; #*** uniqify per index
        }
        else {
            $output_string = "-o $ofile -e $efile";
        }
        
        return qq[$output_string $requirments_string '$cmd$index_spec'];
    }
    
    method determine_queue (VRPipe::Requirements $requirements) {
        #*** supposed to be implemented in eg. VRPipe::Schedulers::LSF if
        #    $self->type eq 'LSF'; hard-coded to something LSF&Sanger specific
        #    for now
        return 'normal';
    }
    
    method run_on_node (Persistent :$submission?, Persistent :$array?, Maybe[PositiveInt] :$index?) {
        $self->throw("submission or array must be supplied") unless $submission || $array;
        my $input_args = '';
        if ($array) {
            $index = $self->get_1based_index($index);
            $submission = VRPipe::PersistentArray->get(id => $array)->member($index); # *** or avoid possible mismatches by directly getting the Submission with the appropriate _hid and _aid?
            $input_args = "array => $array (index $index)";
        }
        else {
            $submission = VRPipe::Submission->get(id => $submission);
            $input_args = "submission => $submission";
        }
        $self->throw("Could not find a submission corresponding to $input_args") unless $submission;
        
        my $job = $submission->job;
        if ($job->pending) {
            $job->run;
        }
        if ($job->finished) {
            $submission->update_status();
        }
    }
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        #*** supposed to be implemented in eg. VRPipe::Schedulers::LSF if
        #    $self->type eq 'LSF'; hard-coded to something LSF&Sanger specific
        #    for now.
        $index ? return $index : return $ENV{foo};
    }

=pod
    # query and manipulate a particular submission (these methods access the
    # scheduler itself, so should normally be avoided; use Submission and Job
    # methods to track state instead):
    if ($sub->scheduled) {
       # $sub->scheduled means that $sub has a ->sid, which is what the following
       # methods extract to work the answers out. It does intellegent caching of
       # data extracted from the scheduler for a given sid, so all these repeated
       # calls may only result in ~1 system call.
       
       if ($sch->pending($sub)) {
           my $current_queue = $sch->queue($sub);
       }
       elsif ($sch->running($sub)) {
           # perhaps you want to clear out a submission that you know from the
           # job itself has already finished, but the scheduler has gotten
           # 'stuck':
           $sch->kill($sub);
           # this tries a normal kill and waits for the ->sid to be cleared, and
           # failing that tries again with a more severe kill/clear if the
           # scheduler supports such a thing. Eventually it will return boolean
           # depending on if the kill was successful (->sid is no longer visible
           # in the list of sids presented by the scheduler).
           
           # or perhaps you've noticed the run time is approaching the limit you
           # set in Requirements, and you want to switch queues before the
           # scheduler kills your job:
           $sch->switch_queues($sub);
           # this method recalcluates the queue given the current Requirments
           # object associated with $sub, and if this is different from ->queue,
           # then it will attempt a queue switch and return 1 if successful, 0 if
           # no switch was necessary, and -1 if the switch failed.
       }
       elsif ($sch->finished($sub)) {
           # these methods fall back on reading the output files on disc if the
           # scheduler has forgotten about the sid in question:
           if ($sch->done($sub)) {
               # yay!
           }
           elsif ($sch->failed($sub)) {
               my $exit_code = $sch->exit_code($sub);
               if ($sch->ran_out_of_memory) { ... }
               elsif ($sch->ran_out_of_time) { ... }
               elsif ($sch->killed) {
                   # Submission finished due to something like ->kill being called;
                   # ie. the user killed the job deliberly, so a resubmission could
                   # work
               }
               else {
                   # Submission ended due to the Job crashing for some other
                   # unknown reason. You should make sure $sub->update_status() got
                   # called, and then use Submission methods to investigate the
                   # STDOUT and STDERR files.
               }
           }
       }
    }
=cut
    
    __PACKAGE__->make_persistent();
}

1;