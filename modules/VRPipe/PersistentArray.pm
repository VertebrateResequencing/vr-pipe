=head1 NAME

VRPipe::PersistentArray - store lists of VRPipe::Persistent objects

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

class VRPipe::PersistentArray extends VRPipe::Persistent {
    has 'members' => (is => 'rw',
                      isa => ArrayRefOfPersistent,
                      lazy => 1,
                      builder => '_build_members',
                      predicate => '_members_populated');
    
    __PACKAGE__->make_persistent(has_many => [_members => 'VRPipe::PersistentArrayMember']);
    
    around get (ClassName|Object $self: ArrayRefOfPersistent :$members?, Persistent :$id?) {
        $self->throw("both id and members cannot be supplied to get() at the same time") if $id && $members;
        $self->throw("get() needs id or members option") unless $id || $members;
        
        if ($id) {
            # get by id
            return $self->$orig(id => $id);
        }
        else {
            # create a new row, then use the new id to create new
            # PersistentArrayMember rows for each supplied member
            my $array = $self->$orig();
            
            my $index = 0;
            foreach my $member (@{$members}) {
                VRPipe::PersistentArrayMember->get(persistentarray => $array, class => ref($member), class_id => $member->id, array_index => ++$index);
            }
            
            $array->members($members);
            
            return $array;
        }
    }
    
    method _build_members {
        # get all PersistentArrayMember rows with our id, and instantiate the
        # corresponding Persistent objects
        my @members;
        foreach my $array_member ($self->_members) {
            push(@members, $self->_array_member_to_member($array_member));
        }
        
        return \@members;
    }
    
    method _array_member_to_member (VRPipe::PersistentArrayMember $array_member) {
        my $class = $array_member->class;
        return $class->get(id => $array_member->class_id);
    }
    
    method member (PositiveInt $index) {
        if ($index == 0) {
            $self->throw("The array index supplied to member() is 1-based");
        }
        
        if ($self->_members_populated) {
            my $members = $self->members;
            return $members->[$index - 1];
        }
        else {
            foreach my $array_member ($self->_members) {
                if ($array_member->array_index == $index) {
                    return $self->_array_member_to_member($array_member);
                }
            }
        }
    }
    
    method size {
        return scalar $self->_members;
    }
}

1;