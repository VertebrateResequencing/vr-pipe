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
    __PACKAGE__->make_persistent(has_many => [members => 'VRPipe::PersistentArrayMember']);
    
    around get (ClassName|Object $self: Persistent :$id!) {
        return $self->$orig(id => $id);
    }
    around create (ClassName|Object $self: ArrayRefOfPersistent :$members!) {
        # create a new row, then use the new id to create new
        # PersistentArrayMember rows for each supplied member
        my $array = $self->$orig();
        
        my $index = 0;
        my @pam_args;
        foreach my $member (@{$members}) {
            push(@pam_args, { persistentarray => $array, class => ref($member), class_id => $member->id, array_index => ++$index });
        }
        VRPipe::PersistentArrayMember->bulk_create_or_update(@pam_args);
        
        return $array;
    }
    
    method member (PositiveInt $index) {
        if ($index == 0) {
            $self->throw("The array index supplied to member() is 1-based");
        }
        
        my ($pam) = VRPipe::PersistentArrayMember->search({ persistentarray => $self->id, array_index => $index });
        $pam || $self->throw("PersistentArray ".$self->id." does not have a member with index $index");
        
        my $class = $pam->class;
        return $class->get(id => $pam->class_id);
    }
    
    method size {
        return scalar VRPipe::PersistentArrayMember->search({ persistentarray => $self->id });
    }
}

1;