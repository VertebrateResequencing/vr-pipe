
=head1 NAME

VRPipe::KeyVal - store arbitrary key value pairs

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

class VRPipe::KeyVal extends VRPipe::Persistent {
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
    
    # to avoid issues with database reserved words, the first column is called
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
    
    around get (ClassName|Object $self: Persistent :$id?, Str :$keyval_key?, Str :$key?, Str :$val?) {
        if ($id) {
            return $self->$orig(id => $id);
        }
        elsif ($key) {
            $keyval_key = $key;
        }
        return $self->$orig(keyval_key => $keyval_key, val => $val);
    }
    
    around create (ClassName|Object $self: Persistent :$id?, Str :$keyval_key?, Str :$key?, Str :$val?) {
        if ($id) {
            return $self->$orig(id => $id);
        }
        elsif ($key) {
            $keyval_key = $key;
        }
        return $self->$orig(keyval_key => $keyval_key, val => $val);
    }
}

1;
