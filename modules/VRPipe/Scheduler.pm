
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
    
    method ensure_running (Str :$cmd!, VRPipe::Requirements :$requirements!, VRPipe::Interface::BackEnd :$backend!, PositiveInt :$count?, Bool :$append_runner_option_to_cmd?) {
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
                my $runner     = VRPipe::Runner->create(cmd => $cmd, aid => $aid);
                my $runner_id  = $runner->id;
                my $run_args   = $append_runner_option_to_cmd ? "append => q[--runner $runner_id]" : '';
                my $worker_cmd = VRPipe::Interface::CmdLine->vrpipe_perl_e(qq[VRPipe::Runner->get(id => $runner_id)->run($run_args)], $backend->deployment);
                $backend->log("the worker cmd is [$worker_cmd]");
                
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
                $backend->log("the scheduler cmd is [$scheduler_cmd_line]");
                
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
    
    method output_dir (PersistentObject $for) {
        my $root_dir = $self->output_root;
        
        my $hashing_string = ref($for) . '::' . $for->id;
        my @subdirs        = $self->hashed_dirs($hashing_string);
        
        $hashing_string =~ s/\:/_/g;
        my $output_dir = dir($root_dir, @subdirs, $hashing_string);
        
        $output_dir->mkpath;
        return $output_dir;
    }
    
    __PACKAGE__->make_persistent();
}

1;
