
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
    
    around get (ClassName|Object $self: Persistent :$id?, ArrayRefOfPersistent :$members?, Bool :$any_order?) {
        $self->throw("You cannot supply both id and members") if $id && $members;
        
        if ($id) {
            return $self->$orig(id => $id);
        }
        elsif ($members) {
            my %paids;
            my $index = 0;
            foreach my $member (@$members) {
                ++$index;
                foreach my $paid (VRPipe::PersistentArrayMember->get_column_values('persistentarray', { class => ref($member), class_id => $member->id, $any_order ? () : (array_index => $index) })) {
                    $paids{$paid}++;
                }
            }
            
            my $expected_count = @$members;
            foreach my $paid (sort { $a <=> $b } keys %paids) {
                next unless $paids{$paid} == $expected_count;
                next unless VRPipe::PersistentArrayMember->search({ persistentarray => $paid }) == $expected_count;
                return $self->$orig(id => $paid);
            }
            return $self->create(members => $members);
        }
        else {
            $self->throw("PersistentArray->get requires id or members");
        }
    }
    
    around create (ClassName|Object $self: ArrayRefOfPersistent :$members!) {
        # create a new row, then use the new id to create new
        # PersistentArrayMember rows for each supplied member
        my $array = $self->$orig();
        
        my $index = 0;
        my @pam_args;
        foreach my $member (@$members) {
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
        $pam || $self->throw("PersistentArray " . $self->id . " does not have a member with index $index");
        
        return $pam->instantiate;
    }
    
    method member_instances {
        my @return;
        foreach my $member ($self->members) {
            push(@return, $member->instantiate);
        }
        return @return;
    }
    
    method size {
        return scalar VRPipe::PersistentArrayMember->search({ persistentarray => $self->id });
    }
}

1;
