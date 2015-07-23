
=head1 NAME

VRPipe::Schedulers::sge - interface to SGE

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This class provides L<Sub (Son of)
GridEngine|https://arc.liv.ac.uk/trac/SGE>-specific command lines for use by
L<VRPipe::Scheduler>.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Schedulers::sge with VRPipe::SchedulerMethodsRole {
    use DateTime;
    use DateTime::Format::Natural;
    use DateTime::TimeZone;
    our $local_timezone = DateTime::TimeZone->new(name => 'local');
    our $dt_parser = DateTime::Format::Natural->new;
    
    our $memory_limit_type;
    our %possible_memory_limit_types = (
        h_data       => 1,
        h_rss        => 1,
        h_vmem       => 1,
        mem_free     => 1,
        s_data       => 1,
        s_rss        => 1,
        s_vmem       => 1,
        virtual_free => 1
    );
    our @prefered_memory_limit_types = qw(mem_free s_rss s_data s_vmem virtual_free h_rss h_data h_vmem);
    our $time_limit_type;
    our %possible_time_limit_types = (h_rt => 1, s_rt => 1);
    our %queue_time_limits;
    our $pe_name;
    our $initialized = 0;
    
    method initialize {
        return if $initialized;
        
        # go through all the queue configs to see what time limits they have set
        # on them
        my @queues;
        open(my $sqlfh, 'qconf -sql |') || $self->throw("Could not open a pipe from qconf -sql");
        while (<$sqlfh>) {
            chomp;
            push(@queues, $_) if $_;
        }
        close($sqlfh);
        
        foreach my $queue (@queues) {
            open(my $sqfh, "qconf -sq $queue |") || $self->throw("Could not open a pipe from qconf -sq $queue");
            my $time_limit;
            while (<$sqfh>) {
                my ($type, $value) = split;
                next unless exists $possible_time_limit_types{$type};
                
                if (!$time_limit_type && $value ne 'INFINITY') {
                    $time_limit_type = $type;
                }
                
                $value = 31536000 if $value eq 'INFINITY';
                if (!$time_limit || $value < $time_limit) {
                    $time_limit = $value;
                }
            }
            close($sqfh) || $self->throw("Could not close a pipe from qconf -sq $queue");
            $queue_time_limits{$queue} = $time_limit;
        }
        
        # try and work out how the sge admin has set up their queue memory
        # limits and pick the best limit type, preferring mem_free over
        # anything else since it makes the most sense and behaves sensibly.
        # We also figure out time limits here as well, only caring if one is
        # required.
        #*** we don't handle multiple types being required...
        open(my $scfh, 'qconf -sc |') || $self->throw("Unable to open a pipe from qconf -sc");
        my (%requestable, %consumable);
        while (<$scfh>) {
            my ($name, undef, $type, undef, $requestable, $consumable) = split;
            $requestable || next;
            next if $requestable eq 'NO';
            if ($type eq 'MEMORY' && exists $possible_memory_limit_types{$name}) {
                if ($requestable eq 'FORCED') {
                    $memory_limit_type = $name;
                }
                elsif ($consumable) {
                    $consumable{$name} = 1;
                }
                else {
                    $requestable{$name} = 1;
                }
            }
            elsif ($type eq 'TIME' && $requestable eq 'FORCED' && exists $possible_time_limit_types{$name}) {
                $time_limit_type = $name;
            }
        }
        close($scfh) || $self->throw("Unable to close a pipe from qconf -sc");
        
        unless ($memory_limit_type) {
            if (keys %consumable) {
                foreach my $type (@prefered_memory_limit_types) {
                    if (exists $consumable{$type}) {
                        $memory_limit_type = $type;
                        last;
                    }
                }
            }
            elsif (keys %requestable) {
                foreach my $type (@prefered_memory_limit_types) {
                    if (exists $requestable{$type}) {
                        $memory_limit_type = $type;
                        last;
                    }
                }
            }
            
            $self->throw("Could not determine how to request memory given the output of qconf -sc; has this sge installation been configured correctly?") unless $memory_limit_type;
        }
        
        # look at all the configured parallel environments to find one with
        # suitable settings
        open(my $splfh, 'qconf -spl |') || $self->throw("Unable to open a pipe from qconf -spl");
        while (<$splfh>) {
            chomp;
            my $pe = $_ || next;
            open(my $spfh, "qconf -sp $pe |") || $self->throw("Unable to open a pipe from qconf -sp $pe");
            my $matches = 0;
            while (<$spfh>) {
                if (/^allocation_rule\s+\$pe_slots/ || /^job_is_first_task\s+TRUE/i) {
                    $matches++;
                }
            }
            close($spfh);
            if ($matches == 2) {
                $pe_name = $pe;
                last;
            }
        }
        close($splfh);
        
        $initialized = 1;
    }
    
    method submit_command (VRPipe::Requirements :$requirements!, Str|File :$stdo_file!, Str|File :$stde_file!, Str :$cmd!, PositiveInt :$count = 1, Int :$global_max = 0, Str :$cwd?) {
        # access the requirements object and build up the string based on
        # memory, time & cpu.
        my $megabytes           = $requirements->memory;
        my $requirements_string = "-l $memory_limit_type=${megabytes}M";
        if ($time_limit_type) {
            my $seconds = $requirements->time;
            $requirements_string .= " -l $time_limit_type=$seconds";
        }
        my $cpus = $requirements->cpus;
        if ($cpus > 1 && $pe_name) {
            $requirements_string .= " -pe $pe_name $cpus";
        }
        
        # deal with job arrays
        my $job_array_args = '';
        if ($count > 1) {
            $job_array_args = " -t 1-$count";
        }
        
        # for command_status() to work efficiently we must always set a job
        # name that corresponds to the cmd. It should also be unique in case of
        # issues with multiple different arrays with the same name
        my $job_name = $self->_job_name($cmd, unique => 1);
        
        #*** we can't cope with complex $cmds that include quotes or semicolons or pipes etc...
        return qq[qsub -N $job_name -o $stdo_file -e $stde_file -m n $requirements_string$job_array_args -V -cwd -b yes $cmd];
    }
    
    method determine_queue (VRPipe::Requirements $requirements, Int $global_max = 0) {
        # we rely on SGE's auto-queue selection, so aren't trying to determine
        # the best queue to run in; determine_queue() is called by a handler for
        # the purpose of queue switching, which we don't support anyway, and its
        # also called via queue_time by a handler; either way we can just return
        # the queue we're currently running in
        return $ENV{QUEUE};
    }
    
    method queue_time (VRPipe::Requirements $requirements) {
        my $queue = $self->determine_queue($requirements);
        return $queue_time_limits{$queue} || 31536000;
    }
    
    method switch_queue (PositiveInt $sid, Str $new_queue) {
        $self->throw("queue switching is not supported and should not be necessary under the SGE scheduler");
    }
    
    method get_scheduler_id {
        return $ENV{JOB_ID};
    }
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        return $index if $index;
        my $aid = $ENV{SGE_TASK_ID};
        if (!$aid || $aid !~ /^\d+$/) {
            return 0;
        }
        return $aid;
    }
    
    method kill_sids (ArrayRef $sid_aids) {
        my @sids;
        foreach my $sid_aid (@$sid_aids) {
            my ($sid, $aid) = @$sid_aid;
            my $id = $aid ? "$sid.$aid" : $sid;
            push(@sids, $id);
            
            if (@sids == 500) {
                system("qdel @sids");
                @sids = ();
            }
        }
        
        system("qdel @sids") if @sids;
    }
    
    method sid_status (PositiveInt $sid, Int $aid) {
        # hmmm, something like 'qstat -j $sid.$aid' does not just give
        # explicit run info on the desired aid, but there is a 'usage' line that
        # seems to indicate if its running
        open(my $qfh, "qstat -j $sid |") || $self->warn("Could not open pipe from qstat -j $sid");
        my $status;
        if ($qfh) {
            $aid ||= 1;
            while (<$qfh>) {
                if (/^usage\s+$aid:\s+/) {
                    $status = 'RUN';
                }
            }
            close($qfh); # (when a $sid does not exist, the close fails; no need to warn about it)
        }
        
        return $status || 'UNKNOWN'; # *** needs to return a word in a defined vocabulary suitable for all schedulers
    }
    
    method run_time (PositiveInt $sid, Int $aid) {
        # sadly something like 'qstat -j $sid.$aid' only gives a submission
        # time, showing info for all members of the array. We'll have to parse
        # the entire output of 'qstat -g d', which shows start time of each
        # array member (if it has started running, otherwise it is the
        # submission time)
        #*** performance issue here, since run_time() may get called many times
        # in a row in loop; we'd like to cache all the start times of all jobs
        # perhaps?
        open(my $qfh, "qstat -g d |") || $self->warn("Could not open pipe from qstat -g d");
        my ($start_epoch, $end_epoch);
        if ($qfh) {
            while (<$qfh>) {
                my ($id, undef, undef, undef, $state, $date, $time, undef, undef, $index) = split;
                next unless $id && $id =~ /^\d+$/;
                next unless $id == $sid;
                if ($aid) {
                    next unless $index && $index == $aid;
                }
                next unless $state eq 'r';
                
                my $dt = $dt_parser->parse_datetime("$date $time");
                $start_epoch = $dt->epoch;
                last;
            }
            close($qfh) || $self->warn("Could not close pipe from qstat -g d");
        }
        
        $start_epoch || return 0;
        $end_epoch ||= time();
        return $end_epoch - $start_epoch;
    }
    
    method command_status (Str :$cmd, PositiveInt :$max?) {
        # qstat doesn't show $cmd, so submit_command() gave jobs a job name we
        # can search for. However qstat also doesn't show the untruncated job
        # name, so we parse the -xml output
        my $count            = 0;
        my @running_sid_aids = ();
        my @to_kill;
        my $job_name_prefix = $self->_job_name($cmd);
        open(my $qfh, "qstat -g d -xml |") || $self->warn("Could not open pipe to qstat -g d -xml");
        if ($qfh) {
            my ($matching_job_name, $state, $sid, $aid);
            while (<$qfh>) {
                if (/^\s*<JB_name>$job_name_prefix.+<\/JB_name>/) {
                    $matching_job_name = 1;
                }
                elsif (/^\s*<state>(.+)<\/state>/) {
                    $state = $1;
                }
                elsif (/^\s*<JB_job_number>(\d+)<\/JB_job_number>/) {
                    $sid = $1;
                }
                elsif (/^\s*<tasks>(\d+)<\/tasks>/) {
                    $aid = $1;
                }
                elsif (/^\s*<\/job_list>/) {
                    $matching_job_name || next;
                    next unless $state =~ /r|qw|t|w|s|S/;
                    $count++;
                    
                    $aid ||= 0;
                    my $sid_aid = "$sid\[$aid\]";
                    
                    if ($state eq 'r') {
                        push(@running_sid_aids, $sid_aid);
                    }
                    elsif ($max && $count > $max) {
                        push(@to_kill, $aid ? "$sid.$aid" : $sid);
                    }
                    
                    $matching_job_name = 0;
                    undef $state;
                    undef $sid;
                    undef $aid;
                }
            }
            close($qfh) || $self->warn("Could not open pipe to qstat -g d -xml");
        }
        
        if (@to_kill) {
            system("qdel @to_kill");
        }
        
        return ($count, \@running_sid_aids);
    }
}

1;
