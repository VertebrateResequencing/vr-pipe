
=head1 NAME

VRPipe::SidToSub - holds the correspondence between scheduler id and submission

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A SidToSub is something created on a compute node that notes what id and index
the scheduling system has given that process. It also notes what farm and
requirements id it was created for.

The server then sees this new SidToSub and assigns it a Submission id. The
process on the compute node can then check its SidToSub, find out what
Submission was assigned to it, then run the Submission.

We do this instead of having the process just selecting its own Submission
because that has proven slow and unreliable - the safest most fool-proof way of
doing things is to have the single server process do assignments so that there
is no risk of a Submission being run by more than one process concurrently.

Having the correspondence between sid[aid] and sub is useful because it then
lets us ask the scheduler to kill it if it is no longer needed, or if it seems
to be stuck (or if we just want to manually investigate a currently running sub
by logging into the compute node it is running on).

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

class VRPipe::SidToSub extends VRPipe::Persistent {
    has 'farm' => (
        is     => 'rw',
        isa    => Varchar [64],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1,
    );
    
    has 'req_id' => (
        is     => 'rw',
        isa    => Persistent,
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'sid' => (
        is     => 'rw',
        isa    => IntSQL [8],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'aid' => (
        is     => 'rw',
        isa    => IntSQL [8],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'sub_id' => (
        is          => 'rw',
        isa         => Persistent,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'assignment_time' => (
        is          => 'rw',
        isa         => Datetime,
        coerce      => 1,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    __PACKAGE__->make_persistent();
    
    after sub_id (Maybe[Int] $sub_id?) {
        shift;
        if (@_) {
            $self->assignment_time(DateTime->now());
        }
    }
    
    method time_since_assignment {
        my $assigned = $self->assignment_time || return 0;
        return time() - $assigned->epoch;
    }
}

1;
