
=head1 NAME

VRPipe::PersistentListRole - store an unordered list of VRPipe::Persistent
objects

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

These lists are immutable and insensitive to order.

Both get() and create() will return an existing list if one had previously been
created using the supplied list of members. Both will also create and return a
new list if one had not been previously created.

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

role VRPipe::PersistentListRole {
    use Digest::MD5;
    use VRPipe::Persistent::InMemory;
    
    requires '_member_class';
    requires '_foreign_key';
    requires '_member_key';
    requires 'lookup';
    
    sub _members_to_string {
        my ($self, $members) = @_;
        return 'undef' unless @$members;
        my @ids = sort { $a <=> $b } map { $_->id } @$members;
        return join(',', @ids);
    }
    
    sub _string_to_lookup {
        my ($self, $str) = @_;
        # limit to printable ascii so hexdigest subroutine will work
        $str =~ tr/\x20-\x7f//cd;
        my $dmd5 = Digest::MD5->new();
        $dmd5->add($str);
        return $dmd5->hexdigest;
    }
    
    method _get_list (ClassName|Object $self: $orig, Maybe[Persistent] $id?, Maybe[ArrayRef] $members?) {
        $self->throw("You cannot supply both id and members") if $id && $members;
        
        if ($id) {
            return $self->$orig(id => $id);
        }
        elsif ($members) {
            # have we previously made a list for these $members?
            my $lookup = $self->_string_to_lookup($self->_members_to_string($members));
            
            my $im       = VRPipe::Persistent::InMemory->new();
            my $lock_key = 'PersistentList.' . $lookup;
            $im->block_until_locked($lock_key);
            
            my ($return) = $self->search({ lookup => $lookup }, { rows => 1 });
            if ($return) {
                $im->unlock($lock_key);
                return $return;
            }
            
            my $created = $self->create($self->_member_key . 's' => $members);
            $im->unlock($lock_key);
            return $created;
        }
        else {
            $self->throw("List->get requires id or members");
        }
    }
    
    method _create_list (ClassName|Object $self: $orig, PersistentArrayRef $members!) {
        # create a new row, then use the new id to create new ListMember rows
        # for each supplied member
        my $lookup = $self->_string_to_lookup($self->_members_to_string($members));
        
        my $im       = VRPipe::Persistent::InMemory->new();
        my $lock_key = 'PersistentList.' . $lookup;
        $im->block_until_locked($lock_key);
        
        my $list = $self->$orig(lookup => $lookup);
        
        my @lm_args;
        my $foreign_key = $self->_foreign_key;
        my $member_key  = $self->_member_key;
        foreach my $member (@$members) {
            push(@lm_args, { $foreign_key => $list->id, $member_key => $member->id });
        }
        my $member_class = $self->_member_class;
        $member_class->bulk_create_or_update(@lm_args);
        
        $im->unlock($lock_key);
        return $list;
    }
    
    method _instantiated_members {
        my @return;
        my $foreign_key  = $self->_foreign_key;
        my $member_class = $self->_member_class;
        my $member_key   = $self->_member_key;
        foreach my $member ($member_class->search({ $foreign_key => $self->id }, { prefetch => $member_key })) {
            push(@return, $member->$member_key);
        }
        return @return;
    }
}

1;
