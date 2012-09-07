
=head1 NAME

VRPipe::PersistentArrayMember - an element of a VRPipe::PersistentArray

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

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

class VRPipe::PersistentArrayMember extends VRPipe::Persistent {
    has 'persistentarray' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::PersistentArray'
    );
    
    has 'class' => (
        is     => 'rw',
        isa    => Varchar [64],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'class_id' => (
        is     => 'rw',
        isa    => Persistent,                        # uncoerced, we want this to be a plain int always
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'array_index' => (
        is     => 'rw',
        isa    => IntSQL [8],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    __PACKAGE__->make_persistent();
}

1;
