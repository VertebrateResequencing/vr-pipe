
=head1 NAME

VRPipe::Persistent::Living - a class to extend from if you can be alive or dead

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Persistent::Living extends VRPipe::Persistent {
    use DateTime;
    
    has 'dead_time' => (
        is      => 'rw',
        isa     => PositiveInt,
        default => 60
    );
    
    has 'heartbeat' => (
        is      => 'rw',
        isa     => Datetime,
        coerce  => 1,
        traits  => ['VRPipe::Persistent::Attributes'],
        default => sub { DateTime->now() }
    );
    
    method time_since_heartbeat {
        my $heartbeat = $self->heartbeat;
        my $t         = time();
        return $t - $heartbeat->epoch;
    }
    
    method alive {
        my $elapsed = $self->time_since_heartbeat;
        my $alive   = $elapsed < $self->dead_time;
        unless ($alive) {
            $self->delete;
        }
        return $alive;
    }
    
    method beat_heart {
        my $still_in_db = $self->search({ id => $self->id });
        return 0 unless $still_in_db;
        $self->heartbeat(DateTime->now());
        $self->update;
        return 1;
    }
}

1;
