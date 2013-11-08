
=head1 NAME

VRPipe::Schedulers::lsf - interface to LSF

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This class provides L<LSF
Platform|http://www.platform.com/workload-management/high-performance-computing>-specific
command lines for use by L<VRPipe::Scheduler>.

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

class VRPipe::Schedulers::lsf with VRPipe::SchedulerMethodsRole {
    use DateTime;
    use DateTime::TimeZone;
    use Path::Class;
    use POSIX qw(ceil);
    our $local_timezone = DateTime::TimeZone->new(name => 'local');
    
    our %months = qw(Jan 1
      Feb 2
      Mar 3
      Apr 4
      May 5
      Jun 6
      Jul 7
      Aug 8
      Sep 9
      Oct 10
      Nov 11
      Dec 12);
    our $date_regex = qr/(\w+)\s+(\d+) (\d+):(\d+):(\d+)/;
    our %queues;
    our @sorted_queues;
    our $mem_limit_multiplier;
    
    method initialize {
        # use lsadmin to see what units memlimit (bsub -M) is in
        my $unit_for_limits = `lsadmin showconf lim | grep LSF_UNIT_FOR_LIMITS`;
        if ($unit_for_limits) {
            my ($unit) = $unit_for_limits =~ /=\s*(\w)/;
            if ($unit) {
                if ($unit eq 'M') {
                    $mem_limit_multiplier = 1;
                }
                elsif ($unit eq 'G') {
                    $mem_limit_multiplier = 0.001;
                }
                # elsif 'K' handled below
            }
        }
        # if something unexpected happened, just assume it's KB
        $mem_limit_multiplier ||= 1000;
        
        # parse bqueues -l to figure out what usable queues we have
        open(my $bqlfh, 'bqueues -l |') || $self->throw("Could not open a pipe to bqueues -l");
        my $queue;
        my $next_is_prio        = 0;
        my $looking_at_defaults = 0;
        my $next_is_memlimit    = 0;
        my $next_is_runlimit    = 0;
        my %highest             = (runlimit => 0, memlimit => 0, max => 0, max_user => 0, users => 0, hosts => 0);
        my $update_highest      = sub {
            my ($type, $new_val) = @_;
            return unless $new_val;
            if ($new_val > $highest{$type}) {
                $highest{$type} = $new_val;
            }
        };
        while (<$bqlfh>) {
            if (/^QUEUE: (\S+)/) {
                $queue = $1;
                next;
            }
            next unless $queue;
            
            if (/^PRIO\s+NICE\s+STATUS\s+MAX\s+JL\/U/) {
                $next_is_prio = 1;
            }
            elsif ($next_is_prio) {
                my ($prio, undef, undef, $max, $max_user) = split;
                $queues{$queue}->{prio}     = $prio;
                $queues{$queue}->{max}      = $max unless $max eq '-';
                $queues{$queue}->{max_user} = $max_user unless $max_user eq '-';
                &$update_highest('max',      $queues{$queue}->{max});
                &$update_highest('max_user', $queues{$queue}->{max_user});
                
                $next_is_prio = 0;
            }
            
            elsif (/^DEFAULT LIMITS:/) {
                $looking_at_defaults = 1;
            }
            elsif (/^MAXIMUM LIMITS:|^SCHEDULING PARAMETERS/) {
                $looking_at_defaults = 0;
            }
            elsif (!$looking_at_defaults) {
                if (/MEMLIMIT/) {
                    my $field = 0;
                    foreach my $word (split(/\s+/, $_)) {
                        next unless $word;
                        $field++;
                        last if $word eq 'MEMLIMIT';
                    }
                    $next_is_memlimit = $field;
                }
                elsif ($next_is_memlimit) {
                    my $field = 0;
                    foreach (/(\d+) (\w)/g) {
                        $field++;
                        if ($field == $next_is_memlimit) {
                            my ($val, $unit) = ($1, $2);
                            if ($unit eq 'G') {
                                $val *= 1000;
                            }
                            elsif ($unit eq 'K') {
                                $val /= 1000;
                            }
                            $queues{$queue}->{memlimit} = $val;
                            &$update_highest('memlimit', $val);
                            last;
                        }
                    }
                    $next_is_memlimit = 0;
                }
                elsif (/RUNLIMIT/) {
                    $next_is_runlimit = 1;
                }
                elsif ($next_is_runlimit && /^\s*(\d+)(?:\.0)? min/) {
                    $queues{$queue}->{runlimit} = $1 * 60;
                    &$update_highest('runlimit', $queues{$queue}->{runlimit});
                    $next_is_runlimit = 0;
                }
            }
            
            if (/^(USERS|HOSTS):\s+(.+?)\s*$/) {
                my $type = lc($1);
                my $vals = $2;
                my @vals = split(/\s+/, $vals);
                if ($type eq 'users') {
                    my %users = map { $_ => 1 } @vals;
                    my $me = getlogin || getpwuid($<);
                    unless (exists $users{all} || exists $users{$me}) {
                        delete $queues{$queue};
                        undef $queue;
                    }
                }
                else {
                    $queues{$queue}->{$type} = scalar(@vals) unless $vals eq 'all';
                    &$update_highest($type, $queues{$queue}->{$type});
                }
            }
            
            if (/^CHUNK_JOB_SIZE:\s+(\d+)/) {
                $queues{$queue}->{chunk_size} = $1;
            }
        }
        close($bqlfh) || $self->throw("Could not close a pipe to bqueues -l");
        
        # for each criteria we're going to sort the queues on later, hard-code
        # [weight, sort-order, significant_change, highest_multiplier]
        # we want to avoid chunked queues because that means jobs will run
        # sequentially instead of in parallel. for time and memory, prefer the
        # queue that is more limited, since we suppose they might be less busy
        # or will at least become free sooner
        my %criteria_handling = (hosts => [10, 1, 25, 2], max_user => [6, 1, 10, 10], max => [5, 1, 20, 5], prio => [4, 1, 50], chunk_size => [10000, 0, 1], runlimit => [2, 0, 3600, 12], memlimit => [2, 0, 16000, 10]);
        
        # fill in some default values for the criteria on all the queues
        my %defaults = (runlimit => 31536000, memlimit => 10000000, max => 10000000, max_user => 10000000, users => 10000000, hosts => 10000000, chunk_size => 0);
        while (my ($criterion, $highest) = each %highest) {
            next unless $highest;
            $defaults{$criterion} = $highest + ($criteria_handling{$criterion}->[2] * $criteria_handling{$criterion}->[3]);
        }
        foreach my $queue (keys %queues) {
            my $queue_hash = $queues{$queue};
            while (my ($criterion, $default) = each %defaults) {
                next if defined $queue_hash->{$criterion};
                $queue_hash->{$criterion} = $default;
            }
        }
        
        # sort the queues, those most likely to run jobs sooner coming first
        my %ranking;
        my %punished_for_max_user;
        foreach my $criterion (qw(max_user max hosts prio chunk_size runlimit memlimit)) { # instead of while loop, because max_user must come first
            my $specs = $criteria_handling{$criterion};
            my ($weight, $sort_order, $significant_change) = @$specs;
            
            my $sort = $sort_order == 0 ? sub { $queues{$a}->{$criterion} <=> $queues{$b}->{$criterion} } : sub { $queues{$b}->{$criterion} <=> $queues{$a}->{$criterion} };
            
            my @sorted = sort $sort keys %queues;
            my $prev_val;
            my $rank = 0;
            foreach my $queue (@sorted) {
                my $val = $queues{$queue}->{$criterion};
                if (defined $prev_val) {
                    my $diff = abs($val - $prev_val);
                    if ($diff >= $significant_change) {
                        if ($criterion eq 'runlimit') {
                            # because the variance in runlimit can be so
                            # massive, increase rank sequentially
                            $rank++;
                        }
                        else {
                            $rank += ceil($diff / $significant_change);
                        }
                    }
                }
                my $punishment = $rank * $weight;
                if ($punishment) {
                    if ($criterion eq 'max_user') {
                        $punished_for_max_user{$queue} = 1;
                    }
                    elsif ($criterion eq 'max') {
                        # don't double-punish for queues that have both max_user
                        # and max
                        $punishment = 0;
                    }
                    $ranking{$queue} += $punishment;
                }
                
                $prev_val = $val;
            }
        }
        
        @sorted_queues = sort { $ranking{$a} <=> $ranking{$b} } keys %ranking;
    }
    
    method submit_command (VRPipe::Requirements :$requirements!, Str|File :$stdo_file!, Str|File :$stde_file!, Str :$cmd!, PositiveInt :$count = 1, Str :$cwd?) {
        # access the requirements object and build up the string based on
        # memory, time, cpu etc.
        my $queue = $self->determine_queue($requirements);
        # *** ...
        my $megabytes          = $requirements->memory;
        my $m                  = $megabytes * $mem_limit_multiplier;
        my $requirments_string = "-q $queue -M$m -R 'select[mem>$megabytes] rusage[mem=$megabytes]'";
        my $cpus               = $requirements->cpus;
        if ($cpus > 1) {
            $requirments_string .= " -n$cpus -R 'span[hosts=1]'";
        }
        
        # for command_status() to work efficiently we must always set a job
        # name that corresponds to the cmd. It must also be unique otherwise
        # LSF would not start running jobs with duplicate names until previous
        # ones complete
        my $job_name = $self->_job_name($cmd, unique => 1);
        
        # work out the scheduler output locations and also specify job arrays
        my $array_def = '';
        my $output_string;
        if ($count > 1) {
            $output_string = "-o $stdo_file.\%I -e $stde_file.\%I";
            $job_name .= qq{\[1-$count]};
        }
        else {
            $output_string = "-o $stdo_file -e $stde_file";
        }
        
        return qq[bsub -J "$job_name" $output_string $requirments_string '$cmd'];
    }
    
    method determine_queue (VRPipe::Requirements $requirements) {
        # pick a queue, preferring ones that are more likely to run our job
        # the soonest
        my $seconds   = $requirements->time;
        my $megabytes = $requirements->memory;
        my $chosen_queue;
        foreach my $queue (@sorted_queues) {
            my $mem_limit = $queues{$queue}->{memlimit};
            next if $mem_limit && $mem_limit < $megabytes;
            my $time_limit = $queues{$queue}->{runlimit};
            next if $time_limit && $time_limit < $seconds;
            
            $chosen_queue = $queue;
            last;
        }
        
        return $chosen_queue;
    }
    
    method queue_time (VRPipe::Requirements $requirements) {
        my $queue = $self->determine_queue($requirements);
        return $queues{$queue}->{runlimit};
    }
    
    method switch_queue (PositiveInt $sid, Str $new_queue) {
        $self->debug("will call bswitch $new_queue $sid");
        system("bswitch $new_queue $sid > /dev/null 2> /dev/null");
    }
    
    method get_scheduler_id {
        return $ENV{LSB_JOBID};
    }
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        $index ? return $index : return $ENV{LSB_JOBINDEX};
    }
    
    method kill_sids (ArrayRef $sid_aids) {
        my @sids;
        foreach my $sid_aid (@$sid_aids) {
            my ($sid, $aid) = @$sid_aid;
            my $id = $aid ? qq{"$sid\[$aid\]"} : $sid;
            push(@sids, $id);
            
            if (@sids == 500) {
                system("bkill @sids");
                @sids = ();
            }
        }
        
        system("bkill @sids") if @sids;
    }
    
    method sid_status (PositiveInt $sid, Int $aid) {
        my $id = $aid ? qq{"$sid\[$aid\]"} : $sid; # when aid is 0, it was not a job array
        open(my $bfh, "bjobs $id |") || $self->warn("Could not call bjobs $id");
        my $status;
        if ($bfh) {
            while (<$bfh>) {
                if (/^$sid\s+\S+\s+(\S+)/) {
                    $status = $1;
                }
            }
            close($bfh);
        }
        
        return $status || 'UNKNOWN';               # *** needs to return a word in a defined vocabulary suitable for all schedulers
    }
    
    method run_time (PositiveInt $sid, Int $aid) {
        my $id = $aid ? qq{"$sid\[$aid\]"} : $sid; # when aid is 0, it was not a job array
        
        open(my $bfh, "bjobs -l $id |") || $self->warn("Could not call bjobs -l $id");
        my ($start_epoch, $end_epoch);
        if ($bfh) {
            my $y = DateTime->now->year;           #*** bad things are going to happen every new year?...
            
            while (<$bfh>) {
                if (/^(.+): .*[Ss]tarted on/) {
                    my $ststr = $1;
                    my ($mo, $d, $h, $m, $s) = $ststr =~ /$date_regex/;
                    my $dt = DateTime->new(year => $y, month => $months{$mo}, day => $d, hour => $h, minute => $m, second => $s, time_zone => $local_timezone);
                    $start_epoch = $dt->epoch;
                }
                elsif (/^(.+): (?:Exited|Completed)/) {
                    my ($mo, $d, $h, $m, $s) = $1 =~ /$date_regex/;
                    my $dt = DateTime->new(year => $y, month => $months{$mo}, day => $d, hour => $h, minute => $m, second => $s, time_zone => $local_timezone);
                    $end_epoch = $dt->epoch;
                    last;
                }
            }
            close($bfh) || $self->warn("Could not call bjobs -l $id");
        }
        
        $start_epoch || return 0;
        $end_epoch ||= time();
        return $end_epoch - $start_epoch;
    }
    
    method command_status (Str :$cmd, PositiveInt :$max?) {
        # bjobs -w does not output a column for both array index and the
        # command. The LSF related modules on CPAN either just parse the command
        # line output or don't work. Ideally we'd use the C-API's
        # lsb_readjobinfo call, but we don't want to be troubled by compilation
        # issues and different versions. We REALLY don't want to manually parse
        # the entire output of bjobs -l for all jobs. Instead when submitting
        # we'll have arranged that JOB_NAME be set to
        # vrpipe_[md5sum_of_cmd]_unique, and in the case of a job array it will
        # have [array_index] appended to it. This lets us use a single
        # bjobs -w call to get all that we need. We can't use -J name to limit
        # which jobs bjobs -w reports on, since we may have submitted the $cmd
        # multiple times as multiple different arrays, each with a uniqified
        # job name. It gets uniquified because otherwise none of the jobs in
        # the second array would start until the first array with the same name
        # ended.
        my $count            = 0;
        my @running_sid_aids = ();
        my @to_kill;
        my $job_name_prefix = $self->_job_name($cmd);
        open(my $bfh, "bjobs -w |") || $self->warn("Could not call bjobs -w");
        if ($bfh) {
            while (<$bfh>) {
                if (my ($sid, $status, $job_name) = $_ =~ /^(\d+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+($job_name_prefix\S+)/) {
                    $job_name || next;
                    next if $status eq 'EXIT';
                    $count++;
                    
                    my ($aid) = $job_name =~ /\[(\d+)\]$/;
                    $aid ||= 0;
                    my $sid_aid = "$sid\[$aid\]";
                    
                    if ($status eq 'RUN') {
                        push(@running_sid_aids, $sid_aid);
                    }
                    elsif ($max && $count > $max) {
                        push(@to_kill, $aid ? qq["$sid_aid"] : $sid);
                    }
                }
            }
            close($bfh) || $self->warn("Could not call bjobs -w");
        }
        
        if (@to_kill) {
            system("bkill @to_kill");
        }
        
        return ($count, \@running_sid_aids);
    }
}

1;
