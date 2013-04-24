
=head1 NAME

VRPipe::SchedulerMethodsRole - a role required by Schedulers

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

L<VRPipe::Scheduler> will look to the site-wide configuration of B<VRPipe> to
see what type of job scheduler should be used, then load
C<VRPipe::Schedulers::[type]> in order to do its work. That class must
implement this role, the required methods of which provide Scheduler the
scheduler-specific command lines needed to do its work.

Until more documentation appears here, see the existing local and lsf types for
a clue about what the required methods are supposed to accept and return.
Contact the author if you get stuck.

The command_status() method, given a cmd to check and maximum count, must
return the count of that cmd currently in the scheduler, and an array ref of
scheduler_id[array_index] strings (array_index == 0 if the scheduler does not
support array refs, or it is not an arrayed job) corresponding to the ones that
are currently running. If more than count jobs for the cmd are in the
scheduler, those that have not yet started to run should be removed from the
scheduler.

run_time(), given a scheduler id and array index, should return how long that
job has actually been running, in seconds. If the job finished already, it
should return the total amount of time spent running.

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

role VRPipe::SchedulerMethodsRole {
    use Digest::MD5;
    
    requires 'start_command';
    requires 'stop_command';
    requires 'submit_command';
    requires 'submit_args';
    requires 'determine_queue';
    requires 'queue_time';
    requires 'switch_queue';
    requires 'get_scheduler_id';
    requires 'get_1based_index';
    requires 'get_sid';
    requires 'kill_sids';
    requires 'sid_status';
    requires 'command_status';
    requires 'run_time';
    
    # this is useful for generating job names based on the cmd
    method _job_name (Str $cmd) {
        my $dmd5 = Digest::MD5->new();
        $dmd5->add($cmd);
        return 'vrpipe_' . $dmd5->hexdigest;
    }
    
    method periodic_method {
        return undef;
    }
    
    method on_exit_method {
        return undef;
    }
}

1;
