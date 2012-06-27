=head1 NAME

VRPipe::Scheduler - a generic interface to job scheduling systems

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

In order to manage the execution of command lines across the many nodes of a
compute cluster, some job scheduling system should be in place. Examples include
LSF and Grid Engine. B<VRPipe> submits the work it wants done to a scheduler;
this Scheduler provides a single consistent interface to all the possible
schedulers.

Currently only LSF is supported. Support for other shedulers can be added by
creating a C<VRPipe::Schedulers::[name]> class that implements the
C<VRPipe::SchedulerMethodsRole> role.

For doing work on the local machine when there isn't a cluster available (or
for testing purposes), there is also a 'Local' scheduler supplied with
B<VRPipe>.

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

class VRPipe::Scheduler extends VRPipe::Persistent {
    use VRPipe::Config;
    my $vrp_config = VRPipe::Config->new();
    use VRPipe::Persistent::SchemaBase;
    use VRPipe::SchedulerMethodsFactory;
    
    has 'type' => (is => 'rw',
                   isa => Varchar[64],
                   builder => 'default_type',
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1,
                   allow_key_to_default => 1);
    
    has 'output_root' => (is => 'ro',
                          isa => Dir,
                          coerce => 1,
                          builder => 'default_output_root',
                          lazy => 1);
    
    method default_type (ClassName|Object $self:) {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment.'_scheduler';
        my $type = $vrp_config->$method_name();
        return lc($type);
    }
    
    method default_output_root (ClassName|Object $self:) {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment.'_scheduler_output_root';
        my $root = $vrp_config->$method_name();
        return "$root"; # stringify what could be a VRPipe::Base::Configuration::Env
    }
    
    # VRPipe::Schedulers::[type] classes will provide scheduler-specific
    # methods
    has 'scheduler_instance' => (is => 'ro',
                                 isa => 'Object',
                                 builder => '_instantiate_method_class',
                                 lazy => 1,
                                 handles => 'VRPipe::SchedulerMethodsRole');
    
    method _instantiate_method_class (ClassName|Object $self:) {
        return VRPipe::SchedulerMethodsFactory->create(lc($self->type), {});
    }
    
    method start_scheduler {
        my $cmd = $self->start_command;
        system("$cmd > /dev/null 2> /dev/null");
    }
    method stop_scheduler {
        my $cmd = $self->stop_command;
        system("$cmd > /dev/null 2> /dev/null");
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

           optionally:
           heartbeat_interval => int (seconds between heartbeats)

=cut
    method submit (VRPipe::Submission :$submission?, VRPipe::Requirements :$requirements?, PersistentArray|ArrayRefOfPersistent :$array?, PositiveInt :$heartbeat_interval?) {
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
            
            if (! $requirements) {
                my %req_ids;
                foreach my $sub (@submissions) {
                    $req_ids{$sub->requirements->id} = 1;
                }
                if (keys %req_ids == 1) {
                    $requirements = $submissions[0]->requirements;
                }
            }
        }
        $self->throw("at least one submission and requirements must be supplied") unless @submissions && $requirements;
        
        # generate a command line that will submit to the scheduler
        my $cmd_line = $self->build_command_line(requirements => $requirements,
                                                 for => $for,
                                                 $heartbeat_interval ? (heartbeat_interval => $heartbeat_interval) : ());
        
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
    
    method build_command_line (VRPipe::Requirements :$requirements!, PersistentObject :$for!, PositiveInt :$heartbeat_interval?) {
        # figure out where STDOUT & STDERR of the scheduler should go
        my $output_dir = $self->output_dir($for);
        
        # the command we want each node to execute once it is run by the
        # scheduler
        my $self_id = $self->id;
        my $node_run_args = $for->isa('VRPipe::PersistentArray') ? "index => shift, array => " : "submission => ";
        $node_run_args .= $for->id;
        if ($heartbeat_interval) {
            $node_run_args .= ", heartbeat_interval => $heartbeat_interval";
        }
        my $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
        my $cmd = qq[perl -MVRPipe::Persistent::Schema -e "VRPipe::Persistent::SchemaBase->database_deployment(q[$deployment]); VRPipe::Scheduler->get(id => $self_id)->run_on_node($node_run_args);"];
        
        return join(' ', $self->submit_command, $self->submit_args(requirements => $requirements,
                                                                   stdo_file => $self->scheduler_output_file($output_dir),
                                                                   stde_file => $self->scheduler_error_file($output_dir),
                                                                   cmd => $cmd,
                                                                   $for->isa('VRPipe::PersistentArray') ? (array => $for) : ()));
    }
    
    method output_dir (PersistentObject $for) {
        my $root_dir = $self->output_root;
        
        my $hashing_string = ref($for).'::'.$for->id;
        my @subdirs = $self->hashed_dirs($hashing_string);
        
        $hashing_string =~ s/\:/_/g;
        my $output_dir = dir($root_dir, @subdirs, $hashing_string);
        
        $output_dir->mkpath;
        return $output_dir;
    }
    
    method scheduler_output_file (Dir $output_dir) {
        return file($output_dir, 'scheduler_stdout');
    }
    method scheduler_error_file (Dir $output_dir) {
        return file($output_dir, 'scheduler_stderr');
    }
    
    method run_on_node (Persistent :$submission?, Persistent :$array?, Maybe[PositiveInt] :$index?, PositiveInt :$heartbeat_interval?) {
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
        if ($heartbeat_interval) {
            $job->heartbeat_interval($heartbeat_interval);
        }
        if ($job->pending) {
            if ($job->block_and_skip_if_ok) {
                # Manager->batch_and_submit does something similar, so this
                # should be unecessary:
                
                # we don't actually block because of race condition issues, and
                # because of fail and restart issues. Instead we always only
                # actually call $job->run if this submission is the first
                # submission created for this $job
                my $schema = $self->result_source->schema;
                my $rs = $schema->resultset('Submission')->search({ 'job' => $job->id }, { rows => 1, order_by => { -asc => 'id' } });
                my $first_sub = $rs->next;
                return unless $first_sub->id == $submission->id;
            }
            $job->run(stepstate => $submission->stepstate);
        }
    }
    
    method wait_for_sid (PositiveInt $sid, Int $aid, PositiveInt $secs = 30) {
        my $t = time();
        while (1) {
            last if time() - $t > $secs;
            
            #*** fork for our call to sid_status and kill child if over time
            #    limit?
            my $status = $self->sid_status($sid, $aid);
            
            last if ($status eq 'UNKNOWN' || $status eq 'DONE' || $status eq 'EXIT');
            sleep(1);
        }
        return 1;
    }
    
    __PACKAGE__->make_persistent();
}

1;