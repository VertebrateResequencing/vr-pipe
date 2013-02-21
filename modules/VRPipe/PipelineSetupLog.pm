
=head1 NAME

VRPipe::PipelineSetupLog - log events that occur to a setup

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

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

class VRPipe::PipelineSetupLog extends VRPipe::Persistent {
    use DateTime;
    
    has 'date' => (
        is     => 'rw',
        isa    => Datetime,
        coerce => 1,
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'ps_id' => (
        is     => 'rw',
        isa    => Persistent,                        # not an actual foreign reference since we want to log stuff about setups that got deleted
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'message' => (
        is     => 'rw',
        isa    => Varchar [254],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'stack' => (
        is          => 'rw',
        isa         => 'Str',
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'de_id' => (
        is                   => 'rw',
        isa                  => IntSQL [9],                        # Int, not Persistent, to allow psuedo-NULL of 0
        traits               => ['VRPipe::Persistent::Attributes'],
        is_key               => 1,
        default              => 0,
        allow_key_to_default => 1
    );
    
    has 'ss_id' => (
        is                   => 'rw',
        isa                  => IntSQL [9],
        traits               => ['VRPipe::Persistent::Attributes'],
        is_key               => 1,
        default              => 0,
        allow_key_to_default => 1
    );
    
    has 'sub_id' => (
        is                   => 'rw',
        isa                  => IntSQL [9],
        traits               => ['VRPipe::Persistent::Attributes'],
        is_key               => 1,
        default              => 0,
        allow_key_to_default => 1
    );
    
    has 'job_id' => (
        is                   => 'rw',
        isa                  => IntSQL [9],
        traits               => ['VRPipe::Persistent::Attributes'],
        is_key               => 1,
        default              => 0,
        allow_key_to_default => 1
    );
    
    __PACKAGE__->make_persistent();
}

1;
