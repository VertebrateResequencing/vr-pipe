
=head1 NAME

VRPipe::Runner - execute a command and know if it is running or not

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A Runner is used to run some command. It lets you check if the command is
running, or if at least it ought to run soon (if it was scheduled to run by a
job management system). It's state tracking beyond that is deliberately simple.
It does not remember if the command you want run already completed successfully
before, nor does it track when and if your command completes, in failure or
otherwise.

You're supposed to use it for commands you simply want to keep running. These
commands should also not care about what directory they're run from (ie. any
paths should be absolute).

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

class VRPipe::Runner extends VRPipe::Persistent::Living {
    use AnyEvent::Util qw(fork_call);
    
    has 'cmd' => (
        is     => 'rw',
        isa    => Text,
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'aid' => (
        is                   => 'rw',
        isa                  => IntSQL [8],
        traits               => ['VRPipe::Persistent::Attributes'],
        default              => 0,
        allow_key_to_default => 1,
        is_key               => 1
    );
    
    has 'sid' => (
        is          => 'rw',
        isa         => IntSQL [8],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'scheduled' => (
        is          => 'rw',
        isa         => Datetime,
        coerce      => 1,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    __PACKAGE__->make_persistent();
    
    method running {
        if ($self->heartbeat && $self->alive) {
            return 1;
        }
        return 0;
    }
    
    after sid (Maybe[Int] $sid?) {
        if ($sid) {
            $self->scheduled(DateTime->now());
        }
    }
    
    method time_scheduled {
        return unless $self->sid;
        my $scheduled = $self->scheduled || return;
        return time() - $scheduled->epoch;
    }
    
    method run (Str :$append?) {
        # start our heartbeat so that another process can check if we're
        # running (it will also stop us running if another process murders us)
        $self->start_beating;
        
        # async fork our cmd
        fork_call {
            my $cmd = $self->cmd;
            $cmd .= ' ' . $append if $append;
            $self->disconnect;
            system($cmd);
            return;
        }
        sub {
            EV::unloop;
        };
        
        # start running
        EV::run;
        
        # now that we've done our job, erase our existence
        $self->end_it_all;
    }
}

1;
