
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
    use EV;
    use AnyEvent;
    
    has 'heartbeat' => (
        is      => 'rw',
        isa     => Datetime,
        coerce  => 1,
        traits  => ['VRPipe::Persistent::Attributes'],
        default => sub { DateTime->now() }
    );
    
    has 'heartbeat_interval' => (
        is      => 'rw',
        isa     => PositiveInt,
        lazy    => 1,
        builder => '_build_default_heartbeat_interval'
    );
    
    has 'heartbeat_timer' => (
        is      => 'rw',
        isa     => 'Object',
        clearer => 'destroy_timer'
    );
    
    has 'dead_time' => (
        is      => 'rw',
        isa     => PositiveInt,
        lazy    => 1,
        builder => '_build_default_dead_time'
    );
    
    method _build_default_heartbeat_interval {
        if (VRPipe::Persistent::SchemaBase->database_deployment eq 'testing') {
            return 3;
        }
        else {
            return 60;
        }
    }
    
    method _build_default_dead_time {
        my $interval = $self->heartbeat_interval;
        
        # try and allow for mysql server time being different to our host time,
        # and other related timing vagaries
        my $multiplier;
        if ($interval < 60) {
            $multiplier = 10;
        }
        elsif ($interval < 300) {
            $multiplier = 3;
        }
        else {
            $multiplier = 1;
        }
        
        return $interval * $multiplier;
    }
    
    method BUILD {
        my $t = time();
        my $timer = EV::periodic 0, $self->heartbeat_interval, 0, sub {
            $self->beat_heart;
        };
        $self->heartbeat_timer($timer);
    }
    
    method time_since_heartbeat {
        my $heartbeat = $self->heartbeat;
        my $t         = time();
        return $t - $heartbeat->epoch;
    }
    
    method alive {
        my $alive = $self->time_since_heartbeat <= $self->dead_time ? 1 : 0;
        unless ($alive) {
            $self->destroy_timer;
            $self->delete;
        }
        return $alive;
    }
    
    method beat_heart {
        my $still_in_db = $self->search({ id => $self->id });
        unless ($still_in_db) {
            $self->destroy_timer;
            return 0;
        }
        $self->heartbeat(DateTime->now());
        $self->update;
        return 1;
    }
}

1;
