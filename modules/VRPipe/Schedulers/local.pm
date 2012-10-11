
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
    my $deployment = VRPipe::Persistent::Schema->database_deployment;
    
    my $ls_script = 'vrpipe-local_scheduler';
    if ($deployment eq 'testing') {
        # we might not have vrpipe-local_scheduler in our PATH, and the
        # modules it needs might not be in our PERL5LIB, so allow
        # us to still work if we're running from the git repo root
        # dir (eg. during testing prior to an install). In fact, for
        # testing purposes, we prefer this to some version of the
        # files installed elsewhere.
        my $local_script = file('scripts', 'vrpipe-local_scheduler');
        my $modules_dir = dir('modules');
        if (-x $local_script && -d $modules_dir) {
            my $thisperl = $Config{perlpath};
            if ($^O ne 'VMS') {
                $thisperl .= $Config{_exe} unless $thisperl =~ m/$Config{_exe}$/i;
            }
            $ls_script = "$thisperl -I$modules_dir $local_script";
        }
    }
    $ls_script .= ' --deployment ' . $deployment;
    
    method start_command {
        return "$ls_script start";
    }
    
    method stop_command {
        return "$ls_script stop";
    }
    
    method submit_command {
        return $ls_script;
    }
    
    method submit_args (VRPipe::Requirements :$requirements!, Str|File :$stdo_file!, Str|File :$stde_file!, Str :$cmd!, VRPipe::PersistentArray :$array?) {
        my $array_def = '';
        my $output_string;
        if ($array) {
            $output_string = "--out $stdo_file.\%I --err $stde_file.\%I";
            my $size = $array->size;
            $array_def = "-a $size ";
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
            
            system("$ls_script kill $id");
            
            sleep(1);
        }
        return 1;
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
}

1;
