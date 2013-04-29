
=head1 NAME

VRPipe::Schedulers::local - stateless 'scheduler' for the local machine

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

The 'Local' scheduler runs jobs directly on the local CPU on a 'when possible'
basis. It holds no state and no queue, so it is possible for a Submission with
high CPU or memory requirements to 'pend' for ages (or even forever if the
requirements exceed the capability of the local machine). Thus this scheduler
is primarily intended for testing purposes only, though could be used in
production on a powerful enough machine. It is recommended to start the server
with --max_submissions set to 1 less than the number of cores in your machine.

It assumes a standard linux type OS with command-line tools such as 'ps' and
'grep', 'kill' and the /proc filesystem. To get this to work on other operating
systems, make a subclass of this and override the following methods as needed:
_build_find_handler_processes_cmd(), _parse_processes(), _kill_cmd(),
_get_process_runtime_cmd(), _process_runtime_to_seconds(),
_build_get_available_cpus_cmd(), _available_cpus_to_count(),
_get_process_memory_cmd(), _process_memory_to_mb(),
_build_get_available_memory_cmd() and _available_memory_to_mb(). If you OS is
very unusual you might also have to override _run_cmd(). Conceptually, however,
this entire module is predicated on the OS being able to provide an id that
uniquely identifies all the processes associated with a running vrpipe-handler.
If it can't do that, a new class may make more sense than a subclass.

You can also subclass this scheduler for use on a cluster of machines. Just
override _cluster_ips to return a list of ips that can be ssh'd to from the
machine the vrpipe-server will run on. It must be possible to do that without
entering a password. Each machine should have the same environment as the
server machine (ie. same environment variables when logged in, same software,
same shared filesystem mounted etc.). It will be blind to anything else running
on those other machines, so do not do this if the other machines are used for
any other purpose or by anyone else. Instead use or write a scheduler module
that interfaces with a 3rd party scheduling system, like LSF.

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

class VRPipe::Schedulers::local with VRPipe::SchedulerMethodsRole {
    use VRPipe::Interface::CmdLine;
    use VRPipe::Persistent::SchemaBase;
    use POSIX qw(ceil);
    
    our $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
    our $backend = VRPipe::Interface::BackEnd->new(deployment => $deployment);
    
    has '_find_handler_processes_cmd' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_build_find_handler_processes_cmd',
    );
    
    has '_get_available_cpus_cmd' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_build_get_available_cpus_cmd',
    );
    
    has '_get_available_memory_cmd' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_build_get_available_memory_cmd',
    );
    
    sub _build_find_handler_processes_cmd {
        return 'ps xo pid,pgid,command | grep vrpipe-handler';
    }
    
    sub _build_get_available_cpus_cmd {
        return qq{perl -MSys::CPU -e 'print Sys::CPU::cpu_count, qq[\n];'};
    }
    
    sub _build_get_available_memory_cmd {
        return 'free -m';
    }
    
    method start_command {
        return 'sleep 1'; #*** not really applicable
    }
    
    method stop_command {
        return 'sleep 1'; #*** not really applicable
    }
    
    method submit_command {
        # we call a method in this module to submit
        my $class = ref($self);
        my ($type) = $class =~ /VRPipe::Schedulers::(\S+)/;
        return VRPipe::Interface::CmdLine->vrpipe_perl_e("VRPipe::Scheduler->get(type => q[$type])->scheduler_instance->submit(\@ARGV)", $deployment);
    }
    
    method submit_args (VRPipe::Requirements :$requirements!, Str|File :$stdo_file!, Str|File :$stde_file!, Str :$cmd!, PositiveInt :$count = 1, Str :$cwd?) {
        # access the requirements object and build up the string based on
        # memory and cpu (other reqs do not apply)
        my $queue              = $self->determine_queue($requirements);
        my $megabytes          = $requirements->memory;
        my $requirments_string = "queue $queue memory $megabytes";
        my $cpus               = $requirements->cpus;
        if ($cpus > 1) {
            $requirments_string .= " cpus $cpus";
        }
        if ($cwd) {
            $requirments_string .= " cwd $cwd";
        }
        
        return qq[$requirments_string count $count cmd '$cmd'];
    }
    
    sub submit {
        my ($self, %args) = @_;
        my $queue     = $args{queue}  || $self->throw("No queue supplied");
        my $megabytes = $args{memory} || $self->throw("No memory supplied");
        my $count     = $args{count}  || $self->throw("No count supplied");
        my $cmd       = $args{cmd}    || $self->throw("No cmd supplied");
        my $cpus      = $args{cpus}   || 1;
        my $cwd       = $args{cwd};
        
        # we do not maintain state or a queue or anything like that. We just
        # check if it is possible to run the $cmd, and do so up to $count times.
        # If we fail to run the command $count times, we don't care; vrpipe-
        # server will just resubmit in the future
        
        # do a pass to see if there's room to run the cmd on any of our ips
        my $usable_ips = $self->_usable_ips($queue, $megabytes, $cpus, $count) || return;
        
        # start a handler on each usable ip (_usable_ips would have returned the
        # same ip more than once if that up could support running more than 1
        # handler)
        $self->_start_new_handlers($usable_ips, $cmd, $cwd ? ($cwd) : ());
    }
    
    method _usable_ips (Str $queue, PositiveInt $megabytes, PositiveInt $cpus, PositiveInt $count) {
        my @usable_ips;
        foreach my $possible_ip ($self->_cluster_ips($queue)) {
            my $available_cpus   = $self->available_cpus($possible_ip);
            my $available_memory = $self->available_memory($possible_ip);
            
            # get the requirements id of every handler process currently running
            my %pgids;
            foreach my $process ($self->_handler_processes($possible_ip)) {
                my ($r) = $process->[2] =~ /-r (\d+) /;
                $pgids{ $process->[1] } = $r || 0; # the setups handler has no -r, and we're going to pretend it uses no cpu, so we can run at least 1 submission handler on single core machines
            }
            
            # get the memory used per handler process group, and don't double-
            # count memory used later on when considering the total available
            # memory
            my %pgid_memory;
            while (my ($pgid, $rid) = each %pgids) {
                # get the total memory used by all processes in this process
                # group
                foreach my $process ($self->_handler_processes($possible_ip, pgid => $pgid)) {
                    my $pid = $process->[0];
                    my $process_memory = $self->memory_used_by_process($possible_ip, $pid);
                    
                    $available_memory -= $process_memory;
                    $pgid_memory{$pgid} += $process_memory;
                }
            }
            
            # calculate how much cpu/memory is in use/reserved
            my $cpus_used   = 0;
            my $memory_used = 0;
            while (my ($pgid, $rid) = each %pgids) {
                my $req = VRPipe::Requirements->get(id => $rid) if $rid;
                $cpus_used += $req ? $req->cpus : 0;
                last if $cpus_used >= $available_cpus;
                
                my $this_memory_used = $pgid_memory{$pgid};
                my $this_memory_reserved = $req ? $req->memory : $this_memory_used;
                if ($this_memory_used > $this_memory_reserved) {
                    # take this opportunity to kill handlers that are using too
                    # much memory
                    $self->_run_cmd($possible_ip, $self->_kill_cmd($pgid));
                    
                    $memory_used += $this_memory_used;
                }
                else {
                    $memory_used += $this_memory_reserved;
                }
                last if $memory_used >= $available_memory;
            }
            
            # see how many times we can submit to this ip
            if ($cpus <= $available_cpus - $cpus_used && $megabytes <= $available_memory - $memory_used) {
                push(@usable_ips, $possible_ip);
                
                while (@usable_ips < $count) {
                    $cpus_used   += $cpus;
                    $memory_used += $megabytes;
                    if ($cpus <= $available_cpus - $cpus_used && $megabytes <= $available_memory - $memory_used) {
                        push(@usable_ips, $possible_ip);
                    }
                    else {
                        last;
                    }
                }
                
                last if @usable_ips == $count;
            }
        }
        
        return \@usable_ips;
    }
    
    method available_cpus (Str $ip) {
        my $available_cpus_str = $self->_run_cmd($ip, $self->_get_available_cpus_cmd) || return 1;
        return $self->_available_cpus_to_count($available_cpus_str);
    }
    
    method _available_cpus_to_count (Str $available_cpus_str) {
        chomp($available_cpus_str);
        return $available_cpus_str;
    }
    
    method available_memory (Str $ip) {
        my $available_memory_str = $self->_run_cmd($ip, $self->_get_available_memory_cmd) || return 0;
        return $self->_available_memory_to_mb($available_memory_str);
    }
    
    method _available_memory_to_mb (Str $free_output) {
        my ($free) = $free_output =~ /cache:\s+\d+\s+(\d+)/;
        if (defined $free) {
            return $free;
        }
        else {
            $backend->log("Unable to get free memory on the system given: $free_output");
            return 0;
        }
    }
    
    method memory_used_by_process (Str $ip, PositiveInt $pid) {
        my $process_memory_str = $self->_run_cmd($ip, $self->_get_process_memory_cmd($pid)) || return 0;
        return $self->_process_memory_to_mb($process_memory_str);
    }
    
    method _get_process_memory_cmd (PositiveInt $pid) {
        return qq[grep VmRSS /proc/$pid/status];
    }
    
    method _process_memory_to_mb (Str $proc_grep) {
        if ($proc_grep =~ /(\d+) kB/) {
            my $grep_bytes = $1 * 1024;
            return ceil($grep_bytes / 1048576);
        }
        else {
            $backend->log("Unable to get process memory given: $proc_grep");
            return 0;
        }
    }
    
    method _start_new_handlers (ArrayRef $usable_ips, Str $cmd, Str $cwd?) {
        my $loop = 0;
        foreach my $ip (@$usable_ips) {
            $loop++;
            my $pgid = $self->_start_new_handler($ip, $cmd, $cwd ? ($cwd) : ());
            if ($pgid) {
                print "Job <$ip:$pgid> is submitted\n";
            }
            else {
                $backend->log("Failed to launch cmd on $ip via ssh");
            }
        }
    }
    
    method _start_new_handler (Str $ip, Str $cmd, Str $working_dir?) {
        # first get the pgids of all handler processes currently running
        my %existing_pgids;
        foreach my $process ($self->_handler_processes($ip)) {
            $existing_pgids{ $process->[1] } = 1;
        }
        
        # now start the new handler process running in the bg
        $self->_run_cmd($ip, $cmd, $working_dir ? (working_dir => $working_dir) : ());
        
        # now work out what the pgid of the new handler must be
        my $pgid;
        foreach my $process ($self->_handler_processes($ip, cmd => $cmd)) {
            my $this_pgid = $process->[1];
            unless (exists $existing_pgids{$this_pgid}) {
                $pgid = $this_pgid;
                last;
            }
        }
        
        return $pgid;
    }
    
    method determine_queue (VRPipe::Requirements $requirements) {
        return 'local';
    }
    
    method queue_time (VRPipe::Requirements $requirements) {
        # we can run for an unlimited time
        return 31536000;
    }
    
    method switch_queue (PositiveInt $sid, Str $new_queue) {
        # we don't support queue switching
        $self->throw("Queue Switching is not supported (and should not be necessary)");
    }
    
    method get_scheduler_id {
        my $pgid = getpgrp(0);
        return "localhost:$pgid";
    }
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        # we don't have any concept of a job 'array', so don't deal with indexes
        return 0;
    }
    
    method kill_sids (ArrayRef $sid_aids) {
        foreach my $sid_aid (@$sid_aids) {
            my ($sid) = @$sid_aid;
            my ($ip, $pgid) = $self->_sid_to_ip_pgid($sid);
            $self->_run_cmd($ip, $self->_kill_cmd($pgid));
        }
    }
    
    method _kill_cmd (PositiveInt $pgid) {
        return qq[kill -TERM -$pgid];
    }
    
    method sid_status (Str $sid, Int $aid) {
        my $processes = $self->_handler_processes($sid, match_pgid => 1);
        return $processes ? 'RUN' : 'UNKNOWN';
    }
    
    method _handler_processes (Str $sid, PositiveInt :$pid?, Bool :$match_pgid?, PositiveInt :$pgid?, Str :$cmd?) {
        my ($ip, $this_pgid) = $self->_sid_to_ip_pgid($sid);
        
        my $processes = $self->_run_cmd($ip, $self->_find_handler_processes_cmd) || return ();
        
        return $self->_parse_processes(
            $processes,
            $pid ? (pid => $pid) : (),
            $match_pgid && $this_pgid ? (pgid => $this_pgid) : (),
            $pgid ? (pgid => $pgid) : (),
            $cmd  ? (cmd  => $cmd)  : ()
        );
    }
    
    method _sid_to_ip_pgid (Str $sid) {
        my ($ip, $pgid) = split(':', $sid);
        $ip ||= $sid;
        return ($ip, $pgid);
    }
    
    method _run_cmd (Str $ip, Str $cmd, Str :$working_dir?) {
        if ($ip eq 'localhost' || $ip eq '127.0.0.1') {
            $cmd = "cd $working_dir; " . $cmd if $working_dir;
            
            if (defined wantarray()) {
                return `$cmd`;
            }
            else {
                # just calling system will not give $cmd its own process group,
                # but we need that; double fork and call setpgrp to be sure
                my $pid = fork;
                if (defined $pid && $pid == 0) {
                    my $pid2 = fork;
                    exit 0 if $pid2;
                    exit 1 if not defined $pid2;
                    setpgrp;
                    $cmd .= ' &>/dev/null &';
                    system($cmd);
                    exit 0;
                }
                waitpid($pid, 0);
            }
        }
        else {
            return $backend->ssh($ip, $cmd, $working_dir ? (working_dir => $working_dir) : ());
        }
    }
    
    method _parse_processes (Str $processes, PositiveInt :$pid?, PositiveInt :$pgid?, Str :$cmd?) {
        my @processes = ();
        foreach my $process (split("\n", $processes)) {
            my ($this_pid, $this_pgid, $this_cmd) = $process =~ /\s*(\d+)\s+(\d+)\s+(.+)$/;
            next if $this_cmd =~ /scheduler_instance->submit/;
            next if $this_cmd =~ /grep vrpipe-handler$/;
            
            if ($pid) {
                next unless $this_pid == $pid;
            }
            if ($pgid) {
                next unless $this_pgid == $pgid;
            }
            if ($cmd) {
                next unless $this_cmd =~ /$cmd/;
            }
            
            push(@processes, [$this_pid, $this_pgid, $this_cmd]);
        }
        
        return @processes;
    }
    
    method run_time (Str $sid, Int $aid) {
        my @processes = $self->_handler_processes($sid, match_pgid => 1);
        my $pid = @processes ? $processes[0]->[0] : return 0;
        
        my ($ip) = $self->_sid_to_ip_pgid($sid);
        my $runtime = $self->_run_cmd($ip, $self->_get_process_runtime_cmd($pid)) || return 0;
        return $self->_process_runtime_to_seconds($runtime);
    }
    
    method _get_process_runtime_cmd (PositiveInt $pid) {
        return qq[ps -p $pid -o etime=];
    }
    
    method _process_runtime_to_seconds (Str $runtime) {
        # [[dd-]hh:]mm:ss
        my ($d, $h, $m, $s) = $runtime =~ /(?:(?:(\d+)-)?(\d+):)?(\d+):(\d+)/;
        $d ||= 0;
        $h ||= 0;
        return ($d * 24 * 60 * 60) + ($h * 60 * 60) + ($m * 60) + $s;
    }
    
    method command_status (Str :$cmd, PositiveInt :$max?) {
        # we have no concept of a pending job, since our submit() method just
        # immediately starts running cmds: we don't have to care about $max or
        # killing things here
        my $count            = 0;
        my @running_sid_aids = ();
        foreach my $ip ($self->_cluster_ips) {
            my %pgids;
            foreach my $process ($self->_handler_processes($ip, cmd => $cmd)) {
                $pgids{ $process->[1] } = 1;
            }
            
            $count += keys %pgids;
            foreach my $pgid (keys %pgids) {
                push(@running_sid_aids, "$ip:$pgid\[0]");
            }
        }
        
        return ($count, \@running_sid_aids);
    }
    
    method _cluster_ips (Str $queue?) {
        return ('localhost');
    }
}

1;
