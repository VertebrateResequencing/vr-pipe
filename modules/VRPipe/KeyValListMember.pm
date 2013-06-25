
=head1 NAME

VRPipe::KeyValListMember - an element of a VRPipe::KeyValList

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

For critical speed reasons, instead of having a 'keyval' column that references
some keyval table, we directly have key and val columns here, so no additional
table lookup is required.

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

class VRPipe::KeyValListMember extends VRPipe::Persistent {
    has 'keyvallist' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::KeyValList'
    );
    
    has 'keyval_key' => (
        is     => 'rw',
        isa    => Varchar [255],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1,
    );
    
    has 'val' => (
        is     => 'rw',
        isa    => Text,
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1,
    );
    
    __PACKAGE__->make_persistent();
    
    # to avoid issues with database reserved words, the key column is called
    # 'keyval_key' instead of just 'key'; let's alias to key
    sub key {
        my $self = shift;
        if (@_) {
            my $key = shift;
            if (length($key)) {
                $self->keyval_key($key);
            }
        }
        return $self->keyval_key;
    }
}

1;
