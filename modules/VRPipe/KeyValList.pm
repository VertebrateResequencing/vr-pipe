
=head1 NAME

VRPipe::KeyValList - store an unordered list of VRPipe::KeyVal objects

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

You probably want to always use get() to safely create or get KeyValLists.

Direct use of create() will always create and return a new KeyValList, even if
one with the same set of keyvals already exists.

These lists are immutable.

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

class VRPipe::KeyValList extends VRPipe::Persistent with VRPipe::PersistentListRole {
    sub _member_class { 'VRPipe::KeyValListMember' }
    sub _member_key   { 'keyval' }
    sub _foreign_key  { 'keyvallist' }
    
    __PACKAGE__->make_persistent(has_many => [members => 'VRPipe::KeyValListMember']);
    
    around get (ClassName|Object $self: Persistent :$id?, ArrayRef[VRPipe::KeyVal] :$keyvals?, HashRef :$hash?) {
        if ($hash) {
            # convert to a list of keyval objects
            $keyvals ||= [];
            my %keyval_ids = map { $_->id => 1 } @$keyvals;
            while (my ($key, $val) = each %{$hash}) {
                my $keyval = VRPipe::KeyVal->create(key => $key, val => $val); # (this returns existing ones if created already)
                next if exists $keyval_ids{ $keyval->id };
                push(@$keyvals, $keyval);
            }
        }
        
        return $self->_get_list($orig, $id, $keyvals);
    }
    
    around create (ClassName|Object $self: ArrayRef[VRPipe::KeyVal] :$keyvals!) {
        return $self->_create_list($orig, $keyvals);
    }
    
    method keyvals {
        return $self->_instantiated_members;
    }
    
    method as_hashref {
        my %return;
        foreach my $keyval ($self->keyvals) {
            $return{ $keyval->key } = $keyval->val;
        }
        return \%return;
    }
}

1;
