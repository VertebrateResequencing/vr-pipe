
=head1 NAME

VRPipe::Schedulers::lsf - interface to LSF

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This class provides L<LSF
Platform|http://www.platform.com/workload-management/high-performance-computing>-specific
command lines for use by L<VRPipe::Scheduler>.

Currently it is unable to start or stop LSF; it is assumed LSF is always up and
running.

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

class VRPipe::Schedulers::lsf with VRPipe::SchedulerMethodsRole {
    our %queues;
    
    method start_command {
        return 'bjobs'; #*** actually, I've no idea how to start lsf
    }
    
    method stop_command {
        return 'bjobs'; #*** actually, I've no idea how to stop lsf
    }
    
    method submit_command {
        return 'bsub';
    }
    
    method submit_args (VRPipe::Requirements :$requirements!, Str|File :$stdo_file!, Str|File :$stde_file!, Str :$cmd!, VRPipe::PersistentArray :$array?) {
        # access the requirments object and build up the string based on memory,
        # time, cpu etc.
        my $queue = $self->determine_queue($requirements);
        # *** ...
        my $megabytes          = $requirements->memory;
        my $m                  = $megabytes * 1000;
        my $requirments_string = "-q $queue -M$m -R 'select[mem>$megabytes] rusage[mem=$megabytes]'";
        my $cpus               = $requirements->cpus;
        if ($cpus > 1) {
            $requirments_string .= " -n$cpus -R 'span[hosts=1]'";
        }
        
        # work out the scheduler output locations and how to pass on the
        # scheduler array index to the perl cmd
        my $index_spec = '';
        my $array_def  = '';
        my $output_string;
        if ($array) {
            $index_spec    = '';                                   #*** something that gives the index to be shifted into perl -e; for LSF we leave it empty and will pick up an env var elsewhere instead
            $output_string = "-o $stdo_file.\%I -e $stde_file.\%I";
            my $size    = $array->size;
            my $uniquer = $array->id;
            $array_def = qq{-J "vrpipe$uniquer\[1-$size]" };
        }
        else {
            $output_string = "-o $stdo_file -e $stde_file";
        }
        
        #*** could solve issue of having to have _aid and _hid in Submission
        #    purely to allow us to find where the lsf output went by having
        #    output go to named pipe that stores directly in db against the
        #    submission id. Same could be done for Job output. Would it work
        #    across the farm?
        
        return qq[$array_def$output_string $requirments_string '$cmd$index_spec'];
    }
    
    method determine_queue (VRPipe::Requirements $requirements) {
        unless (keys %queues) {
            # parse bqueues -l to figure out what usable queues we have
            open(my $bqlfh, 'bqueues -l |') || $self->throw("Could not open a pipe to bqueues -l");
            my $queue;
            my $next_is_prio        = 0;
            my $looking_at_defaults = 0;
            my $next_is_memlimit    = 0;
            my $next_is_runlimit    = 0;
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
                    $queues{$queue}->{max}      = $max eq '-' ? 1000000 : $max;
                    $queues{$queue}->{max_user} = $max_user eq '-' ? 1000000 : $max;
                    
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
                        foreach (/(\d+) K/g) {
                            $field++;
                            if ($field == $next_is_memlimit) {
                                $queues{$queue}->{memlimit} = $1 / 1000;
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
                        $queues{$queue}->{$type} = $vals eq 'all' ? 1000000 : scalar(@vals);
                    }
                }
            }
            close($bqlfh) || $self->throw("Could not close a pipe to bqueues -l");
            
            foreach my $queue (keys %queues) {
                my $queue_hash = $queues{$queue};
                unless (defined $queue_hash->{runlimit}) {
                    $queue_hash->{runlimit} = 31536000;
                }
                unless (defined $queue_hash->{memlimit}) {
                    $queue_hash->{memlimit} = 1000000;
                }
            }
        }
        
        # pick a queue, preferring ones that are more likely to run our job
        # the soonest
        my $seconds   = $requirements->time;
        my $megabytes = $requirements->memory;
        my $chosen_queue;
        foreach my $queue (
            sort {
                     $queues{$b}->{hosts} <=> $queues{$a}->{hosts}
                  || $queues{$b}->{max_user} <=> $queues{$a}->{max_user}
                  || $queues{$b}->{max} <=> $queues{$a}->{max}
                  || $queues{$b}->{prio} <=> $queues{$a}->{prio}
                  || $queues{$a}->{runlimit} <=> $queues{$b}->{runlimit} # for time and memory, prefer the queue that is more limited, since we suppose they might be less busy or will at least become free sooner
                  || $queues{$a}->{memlimit} <=> $queues{$b}->{memlimit}
            } keys %queues
          ) {
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
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        $index ? return $index : return $ENV{LSB_JOBINDEX};
    }
    
    method get_sid (Str $cmd) {
        my $output = `$cmd`;
        my ($sid) = $output =~ /Job \<(\d+)\> is submitted/;
        if ($sid) {
            return $sid;
        }
        else {
            $self->throw("Failed to submit to scheduler");
        }
    }
    
    method kill_sid (PositiveInt $sid, Int $aid, PositiveInt $secs = 30) {
        my $id = $aid ? qq{"$sid\[$aid\]"} : $sid;
        my $t = time();
        while (1) {
            last if time() - $t > $secs;
            
            #*** fork and kill child if over time limit?
            my $status = $self->sid_status($sid, $aid);
            last if ($status eq 'UNKNOWN' || $status eq 'DONE' || $status eq 'EXIT');
            
            system("bkill $id");
            
            sleep(1);
        }
        return 1;
    }
    
    method all_status {
        open(my $bfh, "bjobs |") || $self->warn("Could not call bjobs");
        my %status = ();
        if ($bfh) {
            while (<$bfh>) {
                if (/^(\d+)\s+\S+\s+(\S+)/) {
                    $status{$1} = $2; #*** this does not handle job arrays properly
                }
            }
            close($bfh);
        }
        return %status;
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
}

1;
