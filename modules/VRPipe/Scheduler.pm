
=head1 NAME

VRPipe::Scheduler - a generic interface to job scheduling systems

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

In order to manage the execution of command lines across the many nodes of a
compute cluster, some job scheduling system should be in place. Examples
include LSF and Grid Engine. B<VRPipe> submits the work it wants done to a
scheduler; this Scheduler provides a single consistent interface to all the
possible schedulers.

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
    use VRPipe::Interface::CmdLine;
    use DateTime;
    
    has 'type' => (
        is                   => 'rw',
        isa                  => Varchar [64],
        builder              => 'default_type',
        traits               => ['VRPipe::Persistent::Attributes'],
        is_key               => 1,
        allow_key_to_default => 1
    );
    
    has 'output_root' => (
        is      => 'ro',
        isa     => Dir,
        coerce  => 1,
        builder => 'default_output_root',
        lazy    => 1
    );
    
    method default_type (ClassName|Object $self:) {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment . '_scheduler';
        my $type        = $vrp_config->$method_name();
        return lc($type);
    }
    
    method default_output_root (ClassName|Object $self:) {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment . '_logging_directory';
        my $root        = $vrp_config->$method_name();
        return "$root";       # stringify what could be a VRPipe::Base::Configuration::Env
    }
    
    # VRPipe::Schedulers::[type] classes will provide scheduler-specific
    # methods
    has 'scheduler_instance' => (
        is      => 'ro',
        isa     => 'Object',
        builder => '_instantiate_method_class',
        lazy    => 1,
        handles => 'VRPipe::SchedulerMethodsRole'
    );
    
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
    
    method ensure_running (Str :$cmd!, VRPipe::Requirements :$requirements!, VRPipe::Interface::BackEnd :$backend!, PositiveInt $count?) {
        $count ||= 1;
        
        # see if we're already running this $cmd $count times
        my ($example_runner) = VRPipe::Runner->search({ cmd => $cmd });
        my $running = 0;
        if ($example_runner) {
            my $t             = time();
            my $survival_time = $example_runner->survival_time;
            $running = VRPipe::Runner->search({ cmd => $cmd, heartbeat => { '>' => DateTime->from_epoch(epoch => $t - $survival_time) } });
            warn "got running $running vs count $count\n";
            return if $running >= $count;
            
            # well, are we pending in the scheduler at least?
            my $scheduled = VRPipe::Runner->search({ cmd => $cmd, sid => { '!=' => undef }, scheduled => { '!=' => undef } });
            if ($scheduled) {
                my $not_started = VRPipe::Runner->search({ cmd => $cmd, heartbeat => undef });
                
                # if we've been pending for less than 5mins, just assume that
                # everything is OK and that we are 'running'
                my $five_mins_ago = DateTime->from_epoch(epoch => $t - 300);
                my $short_pends = VRPipe::Runner->search({ cmd => $cmd, heartbeat => undef, scheduled => { '>' => $five_mins_ago } });
                warn "short pends is $short_pends vs not started $not_started\n";
                return if $short_pends >= $not_started;
                $running += $short_pends;
                
                # we think we've been pending in the scheduler for over 5mins
                # now, but are we really?
                my %status = $self->all_status;
                my $long_pend_pager = VRPipe::Runner->search_paged({ cmd => $cmd, heartbeat => undef, scheduled => { '<=' => $five_mins_ago } });
                while (my $runners = $long_pend_pager->next) {
                    foreach my $runner (@$runners) {
                        my $sid = $runner->sid;
                        
                        my $status = $status{$sid} || 'UNKNOWN';
                        if ($status =~ /PEND|RUN/) {
                            if ($status eq 'RUN') {
                                # if the scheduler thinks we're running, give it
                                # 10s and recheck if we're running
                                my $t                = time();
                                my $actually_running = 0;
                                while (1) {
                                    last if time() - $t > 10;
                                    if ($runner->running) {
                                        $actually_running = 1;
                                        $running++;
                                        last;
                                    }
                                    sleep(1);
                                }
                                
                                unless ($actually_running) {
                                    # hmmm, there's something strange going on, just
                                    # kill it
                                    $backend->log("The scheduler told us that $sid was running, but the corresponding Runner had no heartbeat!");
                                    $self->kill_sid($sid, 0, 5);
                                    $runner->delete;
                                }
                            }
                            else {
                                # ok, we're pending, but emit a warning if we've
                                # been pending for over 6hrs in case there's
                                # some kind of problem with the scheduler/ farm
                                #*** email admin only once about this
                                if ($runner->time_scheduled > 21600) {
                                    $backend->log("Runner for cmd [$cmd] has been pending in the scheduler (as $sid) for over 6hrs - is everything OK?");
                                }
                                $running++;
                            }
                        }
                        else {
                            # either we never actually got scheduled, or there's
                            # some problem with the scheduled job, so kill it
                            $self->kill_sid($sid, 0, 5);
                            $runner->delete;
                        }
                    }
                }
            }
            else {
                $backend->log("Runners for cmd [$cmd] were neither running nor scheduled!");
                VRPipe::Runner->search_rs({ cmd => $cmd })->delete;
            }
            
            return if $running >= $count;
        }
        
        my %aids_to_skip;
        if ($running > 0) {
            %aids_to_skip = map { $_ => 1 } VRPipe::Runner->get_column_values('aid', { cmd => $cmd });
            return if keys %aids_to_skip >= $count;
        }
        
        # submit new Runners
        foreach my $aid (0 .. ($count - 1)) {
            next if delete $aids_to_skip{$aid};
            
            my $transaction = sub {
                # the cmd we submit needs a heartbeat, so that when we check if
                # the command is running we can avoid querying the scheduler as
                # much as possible - we wrap $cmd in a Runner
                my $runner = VRPipe::Runner->create(cmd => $cmd, aid => $aid);
                my $worker_cmd = VRPipe::Interface::CmdLine->vrpipe_perl_e(qq[VRPipe::Runner->get(cmd => q{$cmd}, aid => $aid)->run], $backend->deployment);
                
                # construct the command line that will submit our $cmd to the
                # scheduler
                my $scheduler_cmd_line = join(
                    ' ',
                    $self->submit_command,
                    $self->submit_args(
                        requirements => $requirements,
                        stdo_file    => '/dev/null',
                        stde_file    => '/dev/null',
                        cmd          => $worker_cmd
                    )
                );
                
                # go ahead and submit it, getting back an id from the scheduler
                my $sid = $self->get_sid($scheduler_cmd_line);
                
                #*** issues here should also email admin once...
                if ($sid) {
                    # store the sid and scheduled time on the Runner for this
                    # $cmd
                    $runner->sid($sid);
                    $runner->update;
                }
                else {
                    die "The scheduler did not let us submit a job using this command line:\n$scheduler_cmd_line\n";
                }
            };
            $self->do_transaction($transaction, "Failed to submit a Runner");
        }
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
    
    method submit (VRPipe::Submission :$submission?, VRPipe::Requirements :$requirements?, ArrayRefOfPersistent :$array?, PositiveInt :$heartbeat_interval?) {
        my $submissions;
        my $for;
        my $aid = -1;
        if ($submission) {
            $submissions = [$submission];
            $requirements ||= $submission->requirements;
            $for = $submission;
        }
        elsif ($array) {
            my $parray = VRPipe::PersistentArray->create(members => $array);
            $for = $parray;
            $aid++;
            $submissions = $array;
            
            if (!$requirements) {
                my %req_ids;
                foreach my $sub (@$submissions) {
                    $req_ids{ $sub->requirements->id } = 1;
                }
                if (keys %req_ids == 1) {
                    $requirements = $submissions->[0]->requirements;
                }
            }
        }
        $self->throw("at least one submission and requirements must be supplied") unless @$submissions && $requirements;
        
        # generate a command line that will submit to the scheduler
        my $cmd_line = $self->build_command_line(
            requirements => $requirements,
            for          => $for,
            $heartbeat_interval ? (heartbeat_interval => $heartbeat_interval) : ()
        );
        
        # claim all submission objects, associating them with the hashing id,
        # then attempt the submit and set their sid on success or release on
        # failure
        my $for_id      = $for->id;
        my $transaction = sub {
            my $all_claimed = 1;
            foreach my $sub (@$submissions) {
                my $claimed = $sub->claim;
                unless ($claimed) {
                    $all_claimed = 0;
                    last;
                }
                $sub->_hid($for_id);
                $sub->_aid(++$aid);
            }
            
            if ($all_claimed) {
                my $got_sid = $self->get_sid($cmd_line);
                
                if ($got_sid) {
                    foreach my $sub (@$submissions) {
                        $sub->sid($got_sid);
                    }
                    return $got_sid;
                }
                else {
                    foreach my $sub (@$submissions) {
                        $sub->release;
                    }
                    die "failed to submit to scheduler";
                }
            }
            else {
                die "failed to claim all submissions";
            }
        };
        my $sid = $self->do_transaction($transaction, "Failed to claim & submit");
        
        return $sid;
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
        my $cmd        = qq[perl -MVRPipe::Persistent::Schema -e "VRPipe::Persistent::SchemaBase->database_deployment(q[$deployment]); VRPipe::Scheduler->create(id => $self_id)->run_on_node($node_run_args);"];
        
        return join(
            ' ',
            $self->submit_command,
            $self->submit_args(
                requirements => $requirements,
                stdo_file    => $self->scheduler_output_file($output_dir),
                stde_file    => $self->scheduler_error_file($output_dir),
                cmd          => $cmd,
                $for->isa('VRPipe::PersistentArray') ? (array => $for) : ()
            )
        );
    }
    
    method output_dir (PersistentObject $for) {
        my $root_dir = $self->output_root;
        
        my $hashing_string = ref($for) . '::' . $for->id;
        my @subdirs        = $self->hashed_dirs($hashing_string);
        
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
            $index      = $self->get_1based_index($index);
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
