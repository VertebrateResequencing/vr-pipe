
=head1 NAME

VRPipe::Parser::lsf - parse LSF scheduler -o output

=head1 SYNOPSIS
    
    use VRPipe::Parser;
    
    # create object, supplying an lsf -o file
    my $pars = VRPipe::Parser->create('lsf', {file => $file});
    
    # get the array reference that will hold the most recently requested record
    my $parsed_record = $pars->parsed_record();
    
    # loop through the file, getting records; records are parsed in reverse
    # from the end of the file, so you get the most recent information first
    while ($pars->next_record()) {
        my $status = $parsed_record->[1];
        # etc.
    }
    
    # and/or just get a stat of interest from the most recently looked at
    # record (or the most recent (last) record if next_record() had not been
    # used):
    my $status = $pars->status();
    my $mem = $pars->memory();
    my $time = $pars->time();
    # etc.

=head1 DESCRIPTION

The L<LSF
Platform|http://www.platform.com/workload-management/high-performance-computing>
is a job scheduling system. When it is used to schedule and run jobs, it
produces messages when the job completes. This parser is for those report files
(C<-o> files).

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Parser::lsf with VRPipe::ParserRole {
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
    our $date_regex = qr/(\w+)\s+(\d+) (\d+):(\d+):(\d+)\s+(\d+)/;

=head2 parsed_record
 
 Title   : parsed_record
 Usage   : my $parsed_record= $obj->parsed_record()
 Function: Get the data structure that will hold the last parsed record
           requested by next_record()
 Returns : array ref, where the elements are:
           [0]  cmd
           [1]  status
           [2]  memory
           [3]  time
           [4]  cpu_time
           [5]  idle_factor
           [6]  queue
           [7]  sid
           [8]  aid (only if the job was part of a job array)
 Args    : n/a

=cut

=head2 next_record
 
 Title   : next_record
 Usage   : while ($obj->next_record()) { # look in parsed_record }
 Function: Parse the next report from the lsf file, starting with the last and
           working backwards.
 Returns : boolean (false at end of output; check the parsed_record for the
           actual information)
 Args    : n/a

=cut
    
    method next_record {
        # just return if no file set
        my $fh = $self->fh() || return;
        
        # we're only interested in the small LSF-generated report, which may be
        # prefaced by an unlimited amount of output from the program that LSF
        # ran, so we go line-by-line to find our little report. Typically we're
        # only interested in the last result, so we're also actually reading the
        # file backwards
        my ($found_report_start, $found_report_end, $next_is_cmd);
        my ($started, $finished, $cmd, $mem, $status, $queue, $sid, $aid);
        my $cpu = 0;
        while (<$fh>) {
            if (/^Sender: LSF System/) {
                $found_report_start = 1;
                last;
            }
            elsif (/^The output \(if any\) is above this job summary/) {
                $found_report_end = 1;
                next;
            }
            
            if ($found_report_end) {
                if    (/^Started at \S+ (.+)$/)                   { $started     = $1; }
                elsif (/^Subject: Job (\d+)(?:\[(\d+)\])?/)       { $sid         = $1; $aid = $2 if $2; }
                elsif (/^Job was executed.+in queue \<([^>]+)\>/) { $queue       = $1; }
                elsif (/^Results reported at \S+ (.+)$/)          { $finished    = $1; }
                elsif (/^# LSBATCH: User input/)                  { $next_is_cmd = 0; }
                elsif ($next_is_cmd) {
                    $cmd .= $_;
                }
                elsif (!$cmd && $status && /^------------------------------------------------------------/) {
                    $next_is_cmd = 1;
                }
                elsif (/^Successfully completed/)    { $status = 'OK'; }
                elsif (/^Cannot open your job file/) { $status = 'unknown'; }
                elsif (!$status && /^Exited with exit code/) { $status = 'exited'; }
                elsif (/^TERM_\S+ job killed by/)  { $status = 'killed'; }
                elsif (/^TERM_([^:]+):/)           { $status = $1; }
                elsif (/^\s+CPU time\s+:\s*(\S+)/) { $cpu    = $1; }
                elsif (/^\s+Max Memory\s+:\s*(\S+)\s+(\S+)/) {
                    $mem = $1;
                    if ($2 eq 'KB') { $mem /= 1024; }
                    elsif ($2 eq 'GB') { $mem *= 1024; }
                }
            }
        }
        
        # if we didn't see a whole LSF report, assume eof
        unless ($found_report_start && $found_report_end) {
            return;
        }
        unless ($status) {
            $self->warn("a status was not parsed out of a result in " . $self->file);
        }
        
        chomp($cmd) if $cmd;
        
        # calculate wall time and idle factor
        my ($smo, $sd, $sh, $sm, $ss, $sy) = $started  =~ /$date_regex/;
        my ($emo, $ed, $eh, $em, $es, $ey) = $finished =~ /$date_regex/;
        my $dt = DateTime->new(year => $sy, month => $months{$smo}, day => $sd, hour => $sh, minute => $sm, second => $ss);
        my $st = $dt->epoch;
        $dt = DateTime->new(year => $ey, month => $months{$emo}, day => $ed, hour => $eh, minute => $em, second => $es);
        my $et   = $dt->epoch;
        my $wall = $et - $st;
        my $idle = sprintf("%0.2f", ($cpu < 1 ? 1 : $cpu) / ($wall < 1 ? 1 : $wall));
        
        # fill in the parsed_record
        my $pr = $self->parsed_record;
        $pr->[0] = $cmd;
        $pr->[1] = $status;
        $pr->[2] = $mem;
        $pr->[3] = $wall;
        $pr->[4] = $cpu;
        $pr->[5] = $idle;
        $pr->[6] = $queue;
        $pr->[7] = $sid;
        $pr->[8] = $aid;
        
        return 1;
    }

=head2 status
 
 Title   : status
 Usage   : my $status = $obj->status();
 Function: Get the status of the current record, or the last record if
           next_record() hasn't been called yet.
 Returns : string (OK|exited|killed|MEMLIMIT|RUNLIMIT)
 Args    : n/a

=cut
    
    method status {
        return $self->_get_result(1);
    }
    
    method _get_result (Int $index where {$_ >= 0 && $_ <= 8}) {
        my $pr = $self->parsed_record;
        unless (@$pr == 9) {
            $self->next_record || return;
        }
        return $pr->[$index];
    }

=head2 time
 
 Title   : time
 Usage   : my $time = $obj->time();
 Function: Get the wall-time of the current record, or the last record if
           next_record() hasn't been called yet.
 Returns : int (s)
 Args    : n/a

=cut
    
    method time {
        return $self->_get_result(3);
    }

=head2 cpu_time
 
 Title   : cpu_time
 Usage   : my $time = $obj->cpu_time();
 Function: Get the cpu-time of the current record, or the last record if
           next_record() hasn't been called yet.
 Returns : real number (s)
 Args    : n/a

=cut
    
    method cpu_time {
        return $self->_get_result(4);
    }

=head2 idle_factor
 
 Title   : idle_factor
 Usage   : my $idle_factor = $obj->idle_factor();
 Function: Compare cpu time to wall time to see what proportion was spent
           waiting on disc.
 Returns : real number 0-1 inclusive (1 means no time spent waiting on disc, 0
           means the cpu did nothing and we spent all time waiting on disc)
 Args    : n/a

=cut
    
    method idle_factor {
        return $self->_get_result(5);
    }

=head2 memory
 
 Title   : memory
 Usage   : my $memory = $obj->memory();
 Function: Get the max memory used of the current record, or the last record if
           next_record() hasn't been called yet.
 Returns : int (s)
 Args    : n/a

=cut
    
    method memory {
        return $self->_get_result(2);
    }

=head2 cmd
 
 Title   : cmd
 Usage   : my $cmd = $obj->cmd();
 Function: Get the command-line of the current record, or the last record if
           next_record() hasn't been called yet.
 Returns : string
 Args    : n/a

=cut
    
    method cmd {
        return $self->_get_result(0);
    }

=head2 queue
 
 Title   : queue
 Usage   : my $queue = $obj->queue();
 Function: Get the command-line of the current record, or the last record if
           next_record() hasn't been called yet.
 Returns : string
 Args    : n/a

=cut
    
    method queue {
        return $self->_get_result(6);
    }

=head2 sid
 
 Title   : sid
 Usage   : my $sid = $obj->sid();
 Function: Get the submission id of the current record, or the last record if
           next_record() hasn't been called yet.
 Returns : string
 Args    : n/a

=cut
    
    method sid {
        return $self->_get_result(7);
    }

=head2 aid
 
 Title   : aid
 Usage   : my $aid = $obj->aid();
 Function: Get the array index of the current record, or the last record if
           next_record() hasn't been called yet.
 Returns : string
 Args    : n/a

=cut
    
    method aid {
        return $self->_get_result(8);
    }
}

1;
