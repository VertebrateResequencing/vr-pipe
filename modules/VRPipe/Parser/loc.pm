
=head1 NAME

VRPipe::Parser::loc - parse VRPipe local scheduler output

=head1 SYNOPSIS
    
    use VRPipe::Parser;
    
    # create object, supplying a local scheduler file
    my $pars = VRPipe::Parser->create('loc', {file => $file});
    
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
    my $mem = $pars->memory(); # always returns 0, since memory isn't tracked
    my $time = $pars->time();
    # etc.

=head1 DESCRIPTION

When the B<VRPipe> local scheduler is used to schedule and run jobs, it
produces messages when the job completes. This parser is for those report
files.

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

class VRPipe::Parser::loc with VRPipe::ParserRole {
    our $date_regex = qr/(\d+)-(\d+)-(\d+)T(\d+):(\d+):(\d+)/;

=head2 parsed_record
 
 Title   : parsed_record
 Usage   : my $parsed_record= $obj->parsed_record()
 Function: Get the data structure that will hold the last parsed record
           requested by next_record()
 Returns : array ref, where the elements are:
           [0]  cmd
           [1]  status
           [2]  time
           [3]  sid
           [4]  aid
 Args    : n/a

=cut

=head2 next_record
 
 Title   : next_record
 Usage   : while ($obj->next_record()) { # look in parsed_record }
 Function: Parse the next report from the local file, starting with the last and
           working backwards.
 Returns : boolean (false at end of output; check the parsed_record for the
           actual information)
 Args    : n/a

=cut
    
    method next_record {
        # just return if no file set
        my $fh = $self->fh() || return;
        
        # we're only interested in the small generated report, which may be
        # interrupted by an unlimited amount of output from the program that was
        # run, so we go line-by-line to find our little report. Typically we're
        # only interested in the last result, so we're also actually reading the
        # file backwards
        my ($found_report_start, $found_report_end, $next_is_cmd);
        my ($started, $finished, $cmd, $status, $sid, $aid);
        while (<$fh>) {
            #--- vrpipe-local_scheduler report start ---
            # Job: 2[1]
            # Started at: 2011-09-26T10:09:29
            # Running on: uk10k-1-1-01 for user: sb10 in working dir: /nfs/users/nfs_s/sb10/src/git/VertebrateResequencing/vr-pipe
            # Cmd: perl -MVRPipe::Persistent::Schema -e "VRPipe::Persistent::SchemaBase->database_deployment(q[testing]); VRPipe::Scheduler->get(id => 1)->run_on_node(index => shift, array => 2);"
            # STDERR from this Cmd, if any, appears in /lustre/scratch105/vrpipe/testing/d/9/2/8/VRPipe__PersistentArray__2/scheduler_stderr.1
            # STDOUT from this Cmd, if any, appears below
            #---
            #---
            # Finished at 2011-09-26T10:09:49
            # Exit code: 0 (finished OK)
            #--- vrpipe-local_scheduler report end ---
            
            if (/^#--- vrpipe-local_scheduler report start ---/) {
                $found_report_start = 1;
                last;
            }
            elsif (/^#--- vrpipe-local_scheduler report end ---/) {
                $found_report_end = 1;
                next;
            }
            
            if ($found_report_end) {
                if    (/^# Started at: (\S+)$/)       { $started  = $1; }
                elsif (/^# Job: (\d+)(?:\[(\d+)\])?/) { $sid      = $1; $aid = $2 if $2; }
                elsif (/^# Finished at: (\S+)$/)      { $finished = $1; }
                elsif (/^# Cmd: (.+)\n/)              { $cmd      = $1; }
                elsif (/^# Exit code:/) {
                    if    (/OK/)                      { $status = 'OK'; }
                    elsif (/dumped|exited with code/) { $status = 'exited'; }
                    elsif (/killed/)                  { $status = 'killed'; }
                    else                              { $status = 'unknown'; }
                }
            }
        }
        
        # if we didn't see a whole local report, assume eof
        unless ($found_report_start && $found_report_end) {
            return;
        }
        unless ($status) {
            $self->warn("a status was not parsed out of a result in " . $self->file);
        }
        
        chomp($cmd) if $cmd;
        
        # calculate wall time and idle factor
        my ($sy, $smo, $sd, $sh, $sm, $ss) = $started  =~ /$date_regex/;
        my ($ey, $emo, $ed, $eh, $em, $es) = $finished =~ /$date_regex/;
        my $dt = DateTime->new(year => $sy, month => $smo, day => $sd, hour => $sh, minute => $sm, second => $ss);
        my $st = $dt->epoch;
        $dt = DateTime->new(year => $ey, month => $emo, day => $ed, hour => $eh, minute => $em, second => $es);
        my $et   = $dt->epoch;
        my $wall = $et - $st;
        
        # fill in the parsed_record
        my $pr = $self->parsed_record;
        $pr->[0] = $cmd;
        $pr->[1] = $status;
        $pr->[2] = $wall;
        $pr->[3] = $sid;
        $pr->[4] = $aid;
        
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
    
    method _get_result (Int $index where {$_ >= 0 && $_ <= 4}) {
        my $pr = $self->parsed_record;
        unless (@$pr == 5) {
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
        return $self->_get_result(2);
    }
    
    method memory {
        return 0;
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

=head2 sid
 
 Title   : sid
 Usage   : my $sid = $obj->sid();
 Function: Get the submission id of the current record, or the last record if
           next_record() hasn't been called yet.
 Returns : string
 Args    : n/a

=cut
    
    method sid {
        return $self->_get_result(3);
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
        return $self->_get_result(4);
    }
}

1;
