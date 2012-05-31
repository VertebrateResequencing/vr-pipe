=head1 NAME

VRPipe::SchedulerMethodsRole - a role required by Schedulers

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

L<VRPipe::Scheduler> will look to the site-wide configuration of B<VRPipe> to
see what type of job scheduler should be used, then load
C<VRPipe::Schedulers::[type]> in order to do its work. That class must implement
this role, the required methods of which provide Scheduler the
scheduler-specific command lines needed to do its work.

Until more documentation appears here, see the existing local and lsf types
for a clue about what the required methods are supposed to accept and return.
Contact the author if you get stuck.

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

role VRPipe::SchedulerMethodsRole {
    requires 'start_command';
    requires 'stop_command';
    requires 'submit_command';
    requires 'submit_args';
    requires 'determine_queue';
    requires 'switch_queue';
    requires 'get_1based_index';
    requires 'get_sid';
    requires 'kill_sid';
    requires 'sid_status';
}

1;