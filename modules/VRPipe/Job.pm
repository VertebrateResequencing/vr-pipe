
=head1 NAME

VRPipe::Job - state tracking for a command line that must be executed

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A Job describes the actual command line that a pipeline step wishes were
executed. It tracks state related to the execution of its command line, along
with methods to actually run it and get its stdout and stderr.

You do not normally create Job objects yourself directly. A L<VRPipe::Step>
uses one of the C<dispatch*()> methods, and internally this will result in a a
Job being created along with a L<VRPipe::Submission> pointing to it. The
Submission will be submitted to the system's job scheduler, and eventually a
command will run on a node of the compute cluster which gets the Submission,
from which it extracts the Job and calls C<run()> on it.

While a Job executes its command line, it also has another process that emits a
heartbeat, so that B<VRPipe> knows that a Job is running and still healthy. If
something goes wrong on the node and somehow the command fails without the exit
state and fact of completion being recorded normally, the lack of a heartbeat
for a long time will cause B<VRPipe> to consider the Job failed.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2013 Genome Research Limited.

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

class VRPipe::Job extends VRPipe::Persistent {
    use EV;
    use AnyEvent;
    use DateTime;
    use Cwd;
    use Sys::Hostname;
    use Proc::Killfam;
    use VRPipe::Config;
    use Proc::ProcessTable;
    use POSIX qw(ceil);
    use VRPipe::Interface::BackEnd;
    use Scalar::Util qw(weaken isweak);
    
    my $ppt = Proc::ProcessTable->new(cache_ttys => 1);
    our $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
    our $backend = VRPipe::Interface::BackEnd->new(deployment => $deployment);
    
    has 'cmd' => (
        is     => 'rw',
        isa    => Text,
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'dir' => (
        is     => 'rw',
        isa    => Dir,
        coerce => 1,
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'block_and_skip_if_ok' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    has 'exit_code' => (
        is          => 'rw',
        isa         => IntSQL [5],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'pid' => (
        is          => 'rw',
        isa         => IntSQL [6],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'host' => (
        is          => 'rw',
        isa         => Varchar [64],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'user' => (
        is          => 'rw',
        isa         => Varchar [64],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'start_time' => (
        is          => 'rw',
        isa         => Datetime,
        coerce      => 1,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'peak_memory' => (
        is          => 'rw',
        isa         => IntSQL [6],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'end_time' => (
        is          => 'rw',
        isa         => Datetime,
        coerce      => 1,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'output_files' => (
        is      => 'rw',
        isa     => PersistentArrayRef,
        traits  => ['VRPipe::Persistent::Attributes'],
        default => sub { [] }
    );
    
    has 'cmd_monitor' => (
        is      => 'rw',
        isa     => 'Object',
        lazy    => 1,
        builder => '_build_cmd_monitor',
        handles => { start_monitoring => 'start', stop_monitoring => 'stop' }
    );
    
    has '_signalled_to_death' => (
        is  => 'rw',
        isa => 'Str'
    );
    
    has '_watchers' => (
        traits  => ['Array'],
        is      => 'ro',
        isa     => 'ArrayRef[Object]',
        default => sub { [] },
        handles => {
            store_watcher => 'push'
        },
        clearer => 'clear_watchers'
    );
    
    method _build_cmd_monitor {
        my ($watcher, $watcher_count);
        # for the first 2 minutes, when a cmd is most likely to be quickly
        # ramping up its memory usage, we update the peak memory every second
        $watcher = EV::timer_ns 1, 1, sub {
            $self->stop_monitoring if $self->end_time;
            $self->_update_peak_memory;
            
            $watcher_count++;
            if ($watcher_count == 120) {
                # now when memory usage is probably more stable we update every
                # 10 seconds
                $watcher->set(9, 10);
            }
        };
        return $watcher;
    }
    
    method _update_peak_memory (Bool :$no_rounding = 0) {
        my ($host, $pid) = ($self->host, $self->pid);
        unless ($host && $pid) {
            $self->reselect_values_from_db;
            ($host, $pid) = ($self->host, $self->pid);
        }
        
        if ($host && $pid && $host eq hostname()) {
            my $memory = $self->_current_memory($pid);
            if ($memory) {
                # get the stored peak memory; to minimise number of
                # updates we round this up to nearest 100MB
                my $peak_mem = $self->peak_memory;
                if ($peak_mem) {
                    unless ($no_rounding) {
                        my $rounder = 100;
                        my $lower_bound = $peak_mem - ($peak_mem % $rounder);
                        $peak_mem = $lower_bound + $rounder;
                    }
                }
                else {
                    $peak_mem = 0;
                }
                
                if ($memory > $peak_mem) {
                    $self->peak_memory($memory);
                    $self->update;
                }
                if ($memory < 0) {
                    my $stderr_file = $self->stderr_file;
                    my $efh         = $stderr_file->open('>>');
                    print $efh "VRPipe: a bug in Proc::ProcessTable on your system means that we are getting negative RSS values, so we don't really know how much memory this command is using\n";
                    $stderr_file->close;
                }
            }
        }
        $self->disconnect;
    }
    
    method _current_memory (PositiveInt $pid, Bool $include_self = 0) {
        my $pptp;
        foreach my $p (@{ $ppt->table }) {
            if ($p->pid == $pid) {
                $pptp = $p;
                last;
            }
        }
        if ($pptp) {
            my $memory;
            my $pgrp = $pptp->pgrp;
            my %pids_to_skip;
            unless ($include_self) {
                my %parents;
                foreach my $p (@{ $ppt->table }) {
                    if ($p->pgrp == $pgrp) {
                        $parents{ $p->pid } = $p->ppid;
                    }
                }
                
                my $next_to_skip = $$;
                while (exists $parents{$next_to_skip}) {
                    $pids_to_skip{$next_to_skip} = 1;
                    $next_to_skip = $parents{$next_to_skip} || last;
                }
            }
            
            foreach my $p (@{ $ppt->table }) {
                next if exists $pids_to_skip{ $p->pid };
                if ($p->pgrp == $pgrp) {
                    my $rss_bytes = $p->rss;
                    
                    # $p->rss is unreliable, returning negative or too low
                    # values when over about 2GB on some 64bit machines; we'll
                    # have to fall back on assuming this is a typical linux
                    # machine and manually grep /proc
                    my $this_pid = $p->pid;
                    my $grep     = `grep VmRSS /proc/$this_pid/status 2>/dev/null`;
                    my $grep_bytes;
                    if ($grep && $grep =~ /(\d+) kB/) {
                        $grep_bytes = $1 * 1024;
                    }
                    
                    if ($grep_bytes && ($grep_bytes > $rss_bytes)) {
                        $memory += $grep_bytes;
                    }
                    else {
                        $memory += $rss_bytes;
                    }
                }
            }
            
            if ($memory) {
                # return the current memory usage in MB
                return ceil($memory / 1048576);
            }
        }
        return;
    }
    
    __PACKAGE__->make_persistent();
    
    sub _get_args {
        my ($self, $args) = @_;
        unless (exists $args->{dir}) {
            $args->{dir} = cwd();
        }
        else {
            my $dir = dir($args->{dir});
            unless (-d $dir) {
                $dir->mkpath || $self->throw("job output directory '$dir' could not be created");
            }
        }
    }
    
    around get (ClassName|Object $self: %args) {
        $self->_get_args(\%args);
        return $self->$orig(%args);
    }
    
    around create (ClassName|Object $self: %args) {
        $self->_get_args(\%args);
        return $self->$orig(%args);
    }
    
    method ok {
        my $exit_code = $self->exit_code;
        if (defined $exit_code && $exit_code == 0 && $self->end_time) {
            return 1;
        }
        return 0;
    }
    
    method wall_time {
        my $start_time = $self->start_time;
        return 0 unless $start_time;
        my $end_time = $self->end_time;
        if ($end_time) {
            $end_time = $end_time->epoch;
        }
        else {
            $end_time = time();
        }
        return $end_time - $start_time->epoch;
    }
    
    method run (VRPipe::Submission :$submission?, PositiveInt :$allowed_time?) {
        unless ($submission) {
            ($submission) = VRPipe::Submission->search({ job => $self->id, '_done' => 0, '_failed' => 0 }, { rows => 1 });
        }
        our $ss = $submission->stepstate if $submission;
        
        # we'll have 3 responses: -1 = job already exited; 0 = job already
        # running; 1 = we just started running it
        my $response;
        
        # avoid the db-based transaction and locking mechanism below by doing a
        # safer, simple redis lock with NX. We don't rely on this though, since
        # the redis server could go down and we'd lose all locks
        unless ($self->lock) {
            $ss->pipelinesetup->log_event("Job->run() called, but there is a lock on it", dataelement => $ss->dataelement->id, stepstate => $ss->id, submission => $submission->id, job => $self->id) if $ss;
            return 0;
        }
        
        # maintain the lock so that another process can check if we're running
        $self->maintain_lock;
        
        # check we're allowed to run, in a transaction to avoid race condition
        my $start_time  = DateTime->now();
        my $transaction = sub {
            if ($self->start_time) {
                if ($self->end_time) {
                    $response = -1;
                    $ss->pipelinesetup->log_event("Job->run() called, but we already have an end_time", dataelement => $ss->dataelement->id, stepstate => $ss->id, submission => $submission->id, job => $self->id) if $ss;
                }
                else {
                    $ss->pipelinesetup->log_event("Job->run() called, but we've already started running and not finished yet", dataelement => $ss->dataelement->id, stepstate => $ss->id, submission => $submission->id, job => $self->id) if $ss;
                    $response = 0;
                }
                return; # out of the transaction
            }
            
            # set the start time
            $self->reset_job;
            $self->start_time($start_time);
            $self->update;
        };
        $self->do_transaction($transaction, 'Job pending check/ start up phase failed');
        if (defined $response) {
            $self->unlock;
            return $response;
        }
        
        # fork ourselves off a child to run the cmd in. We wrap this up in a
        # single-fire watcher so that nothing actually happens here without the
        # caller doing EV::run
        my $delayed_fork = EV::timer 0, 0, sub {
            weaken($self) unless isweak($self);
            
            my $cmd_pid = fork();
            if (!defined $cmd_pid) {
                $self->unlock;
                $self->throw("attempt to fork cmd failed: $!");
            }
            elsif ($cmd_pid) {
                # parent
            }
            elsif ($cmd_pid == 0) {
                # child, run the command after changing to the correct dir and
                # redirecting output to files
                $self->pid($$);
                $self->host(hostname());
                $self->user(getlogin || getpwuid($<));
                $self->update;
                
                # get all info from db and disconnect before using the info below
                my $dir         = $self->dir;
                my $cmd         = $self->cmd;
                my $stdout_file = $self->stdout_file->path;
                $stdout_file->dir->mkpath;
                my $stderr_file = $self->stderr_file->path;
                $self->disconnect;
                
                open STDOUT, '>', $stdout_file or $self->throw("Can't redirect STDOUT to '$stdout_file': $!");
                untie *STDERR;
                open STDERR, '>', $stderr_file or $self->throw("Can't redirect STDERR to '$stderr_file': $!");
                
                chdir($dir);
                # exec is supposed to get our $cmd to run whilst keeping the
                # same $cmd_pid, but on some systems like Ubuntu the sh (dash)
                # is a bit sucky: http://www.perlmonks.org/?node_id=785284 We
                # can't force list mode in the normal way because we actually
                # require the use of the shell to do things like run multi-line
                # commands and pipes etc. $cmd having a different pid to
                # $cmd_pid matters because we need to know the correct pid if
                # the server needs to kill it later. Instead we force the use of
                # a better shell (default bash) for everything, which might be
                # less efficient in some cases, but the difference is going to
                # be meaningless for us.
                my $shell = VRPipe::Config->new->exec_shell;
                if ($shell) {
                    exec {$shell} $shell, '-c', $cmd;
                }
                else {
                    exec $cmd;
                }
            }
            
            # start our watcher that keeps track of peak memory usage
            $self->start_monitoring;
            
            # setup signal watchers for the various ways a scheduler might try
            # to ask us to die; in response we will try and guess why and update
            # resource requirements on all submissions that use us, then kill
            # the cmd
            foreach my $signal (qw(TERM INT QUIT USR1 USR2)) {
                my $signal_watcher = EV::signal $signal, sub {
                    $self->_update_peak_memory(no_rounding => 1); # the child should actually be killed by now, so this likey does nothing
                    my $own_memory = $self->_current_memory($$, 1);
                    my $total_memory = $own_memory;
                    
                    my $explanation = "VRPipe: the cmd received SIG$signal";
                    
                    # look at the requirements of the submission that was used
                    # to run us
                    my $sid_str = '';
                    if ($submission) {
                        my ($sid_to_sub) = VRPipe::SidToSub->search({ sub_id => $submission->id });
                        if ($sid_to_sub) {
                            my $sid_aid = $sid_to_sub->sid . '[' . $sid_to_sub->aid . ']';
                            $sid_str = "; scheduler id $sid_aid";
                        }
                        
                        my $changed      = 0;
                        my $peak_memory  = $self->peak_memory || 0;
                        my $fudge_factor = 100;
                        $total_memory += $peak_memory + $fudge_factor;
                        my $reserved_memory = $submission->memory;
                        my $mem_cmp         = '<';
                        if ($total_memory >= $reserved_memory) {
                            $mem_cmp = '>';
                            $submission->extra_memory;
                            $changed = 1;
                        }
                        $explanation .= qq[; total memory used (cmd ${peak_memory}MB + vrpipe ${own_memory}MB + ${fudge_factor}MB fudge) $mem_cmp reserved memory (${reserved_memory}MB)];
                        
                        if ($allowed_time) {
                            my $wall_time = $self->wall_time;
                            my $time_cmp  = '<';
                            if ($wall_time >= $allowed_time) {
                                $time_cmp = '>';
                                my $new_time = (ceil($wall_time / 60 / 60) + 1) * 60 * 60;
                                $submission->extra_time($new_time - $submission->time);
                                $changed = 1;
                            }
                            $explanation .= qq[; wall time used (${wall_time}s) $time_cmp time allowed in queue (${allowed_time}s)];
                        }
                        
                        if ($changed) {
                            # also change the requirements for all other subs
                            # for us *** does this happen any more?
                            my $new_req_id = $submission->requirements->id;
                            VRPipe::Submission->search_rs({ job => $self->id, requirements => { '!=' => $new_req_id } })->update({ requirements => $new_req_id });
                        }
                    }
                    
                    my $stderr_file = $self->stderr_file;
                    if ($stderr_file) {
                        my $efh = $stderr_file->open('>>');
                        print $efh $explanation, "\n";
                        $stderr_file->close;
                    }
                    $ss->pipelinesetup->log_event("Job->run() signal watcher detected SIG$signal ($explanation$sid_str), will kill_job()", dataelement => $ss->dataelement->id, stepstate => $ss->id, submission => $submission->id, job => $self->id) if $ss;
                    
                    $self->_signalled_to_death($signal);
                    
                    $self->kill_job($submission);
                    $self->disconnect;
                    $self->unlock;
                };
                $self->store_watcher($signal_watcher);
            }
            
            #*** what we don't handle well is the situation where we get
            #    SIGSTOPed, the node running us becomes unresponsive, this Job
            #    starts up on another node, fails to kill our pid and starts
            #    running cmd again, then this node comes back to life and we are
            #    SIGCONTed. Our loss of lock will cause us to throw, but only
            #    after up to 60 seconds - plenty of time for our cmd_pid to
            #    corrupt the output of the other running Job. I tried
            #    experimenting with a SIGCONT watcher, and with turning trace on
            #    for the child watcher, but these didn't seem to let me kill the
            #    cmd_pid before it resumed
            
            # set up a watcher that detects when the child exits
            my $child_watcher;
            $child_watcher = EV::child $cmd_pid, 0, sub {
                my $exit_code = $child_watcher->rstatus;
                waitpid($cmd_pid, 0); # is this necessary??
                $self->stop_monitoring;
                
                $ss->pipelinesetup->log_event("Job->run() cmd-running child exited with code $exit_code", dataelement => $ss->dataelement->id, stepstate => $ss->id, submission => $submission->id, job => $self->id) if $ss;
                
                if ($exit_code == 0) {
                    # we shouldn't have exited with success yet lost our lock,
                    # so in this case we return early before updating any
                    # columns on job or submission to avoid screwing up another
                    # process that might now be running the same job
                    # (getting killed can lose our lock, and in that case we
                    # will continue and record the submission as failed)
                    if ($submission && !$submission->locked(by_me => 1)) {
                        $self->warn("The Submission's (" . $submission->id . ") lock has been lost prior to completing Job status");
                        return;
                    }
                    if (!$self->locked(by_me => 1)) {
                        $self->warn("Our (job " . $self->id . ") lock has been lost prior to completing Job status");
                        return;
                    }
                }
                
                #*** ideally we want to update_stats_from_disc on all our
                # output files, but that is time consuming and potentially
                # risky: we can end up exiting during that phase before we
                # stop beating and set our end_time and exit_code, which is more
                # critical
                my $end_time = DateTime->now();
                my @to_trigger;
                my $transaction = sub {
                    # update ourselves (our still-maintained redis lock should
                    # ensure we're seeing the latest data)
                    @to_trigger = ();
                    
                    unless ($self->start_time) {
                        #*** we can somehow manage to have no start
                        # time?? This is very bad, but there's not much
                        # we can do about it now other than set it to end_time
                        $self->start_time($end_time);
                    }
                    $self->exit_code($exit_code);
                    $self->end_time($end_time);
                    $self->update;
                    
                    # update our Submission
                    if ($submission) {
                        if ($exit_code == 0) {
                            # say the sub is done
                            $submission->_done(1);
                            $submission->_failed(0);
                            $submission->_claim(0);
                            
                            # we completed successfully; if this was the last
                            # submission of its step and they're all done now,
                            # will do a trigger for this dataelement to create
                            # submissions for the next step; compare counts of
                            # total submissions for our step to count of those
                            # done to see if the step is complete. We can't just
                            # get quick counts though otherwise the for =>
                            # update lock doesn't do anything
                            my $stepstate_id = $ss->id;
                            my @step_subs    = VRPipe::Submission->search({ stepstate => $stepstate_id });
                            my $done         = 0;
                            foreach my $sub (@step_subs) {
                                next if $sub->id == $submission->id;
                                $sub->done || last;
                                $done++;
                            }
                            $submission->update;
                            
                            if ($done == @step_subs - 1) {
                                $ss->pipelinesetup->log_event("At end of Job->run() noted that all Submissions for this Job's Submission's StepState are done, so will trigger the next Step", stepstate => $ss->id, dataelement => $ss->dataelement->id);
                                push(@to_trigger, [$ss->pipelinesetup, $ss->dataelement]);
                                
                                # also trigger any dataelements that have the
                                # same submissions as $step_state
                                my $others = VRPipe::StepState->get_column_values([qw(pipelinesetup dataelement)], { same_submissions_as => $stepstate_id });
                                foreach my $psde (@$others) {
                                    my ($sid, $de) = @$psde;
                                    push(@to_trigger, [VRPipe::PipelineSetup->get(id => $sid), VRPipe::DataElement->get(id => $de)]);
                                }
                            }
                        }
                        else {
                            # say the sub is failed
                            $submission->_done(0);
                            $submission->_failed(1);
                            $submission->_claim(0);
                            $submission->update;
                        }
                    }
                };
                $self->do_transaction($transaction, "Failed to finalise Job and Submission state after the cmd-running child exited");
                
                # trigger the next step now, outside the above transaction,
                # since trigger() isn't really idempotent: it deletes files.
                # Problems with triggering may also explain why we can trigger
                # the next step yet end up not setting job end_time, which can
                # be disastrous if the trigger had a behaviour that deleted our
                # inputs and then the job gets tried again because it looks
                # dead.
                if (@to_trigger) {
                    foreach my $ref (@to_trigger) {
                        my ($setup, $de) = @$ref;
                        my $error_message = $setup->trigger(dataelement => $de);
                        
                        if ($error_message) {
                            $setup->log_event("Job->run() completed and tried to trigger the next step but failed: $error_message", dataelement => $de->id, stepstate => $ss->id, submission => $submission->id, job => $self->id);
                        }
                        else {
                            $setup->log_event("Job->run() completed and successfully triggered the next step", dataelement => $de->id, stepstate => $ss->id, submission => $submission->id, job => $self->id);
                        }
                    }
                }
                
                # update and unlock now (not before the above, since that can
                # take a while, and we don't want another process thinking we
                # died before setting end_time and triggering)
                $self->update;
                $self->disconnect;
                
                #*** theoretically updating file existence now might be too
                # late, but we assume StepRole will recheck file existence
                # on files that seem to be missing anyway
                if ($self->pid) {
                    #*** somehow we can not have a pid, which breaks the
                    # following 2 calls
                    $self->stdout_file->update_stats_from_disc(retries => 3);
                    $self->stderr_file->update_stats_from_disc(retries => 3);
                }
                
                my $o_files = $self->output_files;
                if ($o_files && @$o_files) {
                    foreach my $o_file (@$o_files) {
                        $o_file->update_stats_from_disc(retries => 3);
                    }
                }
                elsif ($submission) {
                    $submission->stepstate->update_output_file_stats;
                }
                
                $self->disconnect;
                $self->unlock;
                $self->clear_watchers;
            };
            $self->store_watcher($child_watcher);
        };
        $self->store_watcher($delayed_fork);
        
        return 1;
    }
    
    method stdout_file {
        my $dir = $self->dir;
        my $file = $self->_std_file || return;
        return VRPipe::File->create(path => file($dir, $file . '.o'), type => 'txt');
    }
    
    method stderr_file {
        my $dir = $self->dir;
        my $file = $self->_std_file || return;
        return VRPipe::File->create(path => file($dir, $file . '.e'), type => 'txt');
    }
    
    method _std_file {
        my $pid = $self->pid || return;
        my $host = $self->host;
        return ".host_$host.pid_$pid";
    }
    
    method kill_job (VRPipe::Submission $submission?) {
        unless ($submission) {
            ($submission) = VRPipe::Submission->search({ job => $self->id }, { rows => 1 });
        }
        my $ss = $submission->stepstate if $submission;
        $ss->pipelinesetup->log_event("Job->kill_job() called", dataelement => $ss->dataelement->id, stepstate => $ss->id, job => $self->id, record_stack => 1) if $ss;
        
        my ($user, $host, $pid) = ($self->user, $self->host, $self->pid);
        if ($user && $host && $pid) {
            if (hostname() eq $host) {
                killfam "KILL", $pid;
            }
            else {
                eval {
                    $backend->ssh("$user\@$host", qq[perl -MProc::Killfam -e 'killfam q[KILL], $pid']);
                    # *** we need to do something if the kill fails...
                };
                if ($@) {
                    die unless $@ eq "ssh timed out\n";
                    # *** and how do we handle not being able to ssh into the host
                    #     at all?
                    $self->warn("ssh to $host timed out, assuming that process $pid is dead...");
                }
            }
            
            my $ofiles = $self->output_files;
            if ($ofiles && @$ofiles) {
                foreach my $file (@$ofiles) {
                    $file->unlink;
                }
            }
            elsif ($submission) {
                # if the Step only called dispatch once (ie. only 1 sub per
                # StepState, Step author may not have specified job output
                # files, leaving it up to automagial handling. See if our sub
                # is the only one for the StepState
                my $ss       = $submission->stepstate;
                my @all_subs = $ss->submissions;
                if (@all_subs == 1) {
                    $ss->unlink_output_files;
                    $ss->unlink_temp_files;
                }
            }
            
            $self->unlock;
        }
    }
    
    method reset_job {
        $self->exit_code(undef);
        $self->pid(undef);
        $self->host(undef);
        $self->user(undef);
        $self->start_time(undef);
        $self->peak_memory(undef);
        $self->end_time(undef);
        $self->update;
        $self->reselect_values_from_db;
        return 1;
    }
}

1;
