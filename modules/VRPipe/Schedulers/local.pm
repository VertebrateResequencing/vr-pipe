
=head1 NAME

VRPipe::Schedulers::local - interface to the Local scheduler

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

The 'Local' scheduler is supplied with B<VRPipe> and runs jobs on the local CPU
with a very simple queuing and scheduling system. It is primarily for testing
purposes only.

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

class VRPipe::Schedulers::local with VRPipe::SchedulerMethodsRole {
    use Config;
    use VRPipe::Persistent::Schema;
    use VRPipe::Interface::CmdLine;
    my $deployment = VRPipe::Persistent::Schema->database_deployment;
    
    my $ls_script = VRPipe::Interface::CmdLine->vrpipe_script_command('vrpipe-local_scheduler', $deployment);
    
    method start_command {
        return "$ls_script start";
    }
    
    method stop_command {
        return "$ls_script stop";
    }
    
    method submit_command {
        return $ls_script;
    }
    
    method submit_args (VRPipe::Requirements :$requirements!, Str|File :$stdo_file!, Str|File :$stde_file!, Str :$cmd!, PositiveInt :$count = 1) {
        my $array_def = '';
        my $output_string;
        if ($count) {
            $output_string = "--out $stdo_file.\%I --err $stde_file.\%I";
            $array_def     = "-a $count ";
        }
        else {
            $output_string = "--out $stdo_file --err $stde_file";
        }
        
        return qq[$array_def$output_string submit '$cmd'];
    }
    
    method determine_queue (VRPipe::Requirements $requirements) {
        return 'local';
    }
    
    method queue_time (VRPipe::Requirements $requirements) {
        return 1314000;
    }
    
    method switch_queue (PositiveInt $sid, Str $new_queue) {
        return;
    }
    
    method get_scheduler_id {
        return $ENV{VRPIPE_LOCAL_JOBID};
    }
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        $index ? return $index : return $ENV{VRPIPE_LOCAL_JOBINDEX};
    }
    
    method get_sid (Str $cmd) {
        my $output = `$cmd`;
        my ($sid) = $output =~ /Job \<(\d+)\> is submitted/;
        if ($sid) {
            return $sid;
        }
        else {
            $self->throw("Failed to submit to scheduler given command $cmd");
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
            
            system("$ls_script kill $id");
            
            sleep(1);
        }
        return 1;
    }
    
    method batch_kill_sids (ArrayRef $sid_aids) {
        # unlike kill_sid(), we're all about speed, so can't care about if the
        # kill actually worked or not
        
        my @sids;
        foreach my $sid_aid (@$sid_aids) {
            my ($sid, $aid) = @$sid_aid;
            my $id = $aid ? qq{"$sid\[$aid\]"} : $sid;
            push(@sids, $id);
            
            if (@sids == 500) {
                system("$ls_script kill @sids");
                @sids = ();
            }
        }
        
        system("$ls_script kill @sids") if @sids;
    }
    
    method all_status {
        open(my $bfh, "$ls_script jobs |") || $self->warn("Could not call $ls_script jobs");
        my %status = ();
        if ($bfh) {
            while (<$bfh>) {
                if (/^(\d+)\s+\S+\s+(\S+)/) {
                    $status{$1} = $2;
                }
            }
            close($bfh);
        }
        return %status;
    }
    
    method sid_status (PositiveInt $sid, Int $aid) {
        my $id = $aid ? qq{"$sid\[$aid\]"} : $sid; # when aid is 0, it was not a job array
        open(my $bfh, "$ls_script jobs $id |") || $self->warn("Could not call vrpipe-local_scheduler jobs $id");
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
        
        open(my $bfh, "$ls_script jobs $id |") || $self->warn("Could not call $ls_script jobs $id");
        my $run_time = 0;
        if ($bfh) {
            while (<$bfh>) {
                if (/^\d+\t\d+\t\S+\t\S+\t\S+\t\S+\t(\d+)/) {
                    $run_time = $1;
                    last;
                }
            }
            close($bfh) || $self->warn("Could not call $ls_script jobs $id");
        }
        
        return $run_time;
    }
    
    method command_status (Str :$cmd, PositiveInt :$max?) {
        my $count            = 0;
        my @running_sid_aids = ();
        my @to_kill;
        open(my $bfh, "$ls_script jobs |") || $self->warn("Could not call $ls_script jobs");
        if ($bfh) {
            while (<$bfh>) {
                if (my ($sid, $aid, $status) = $_ =~ /^(\d+)\t(\d+)\t(\S+)\t\S+\t\S*\t\S+\t\d+\t$cmd/) {
                    $status || next;
                    next if $status eq 'EXIT';
                    $count++;
                    
                    my $sid_aid = "$sid\[$aid\]";
                    if ($status eq 'RUN') {
                        push(@running_sid_aids, $sid_aid);
                    }
                    elsif ($max && $count > $max) {
                        push(@to_kill, $sid_aid);
                    }
                }
            }
            close($bfh) || $self->warn("Could not call $ls_script jobs");
        }
        
        if (@to_kill) {
            system("$ls_script kill @to_kill");
        }
        
        return ($count, \@running_sid_aids);
    }
}

1;
