
=head1 NAME

VRPipe::PersistentListRole - store an unordered list of VRPipe::Persistent
objects

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

You probably want to always use get() to safely create or get Lists.

Direct use of create() will always create and return a new List, even if one
with the same set of members already exists.

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

role VRPipe::PersistentListRole {
    requires '_member_class';
    requires '_foreign_key';
    requires '_member_key';
    
    our %empty_lists;
    
    sub _class {
        my $self = shift;
        return ref($self) ? ref($self) : $self;
    }
    
    method _empty_list (ClassName|Object $self: $orig) {
        my $class = $self->_class;
        
        if (defined $empty_lists{$class}) {
            return $empty_lists{$class};
        }
        
        #*** I have no idea how to do this in DBIx::Class...
        # SELECT *
        # FROM keyvallist
        # LEFT OUTER JOIN keyvallistmember
        #   ON (keyvallistmember.keyvallist = keyvallist.id)
        #   WHERE keyvallistmember.keyvallist IS NULL
        my $member_class = $self->_member_class;
        my $pager = $class->get_column_values_paged('id', {});
        my $empty;
        PAGES: while (my $ids = $pager->next) {
            foreach my $l_id (@$ids) {
                my $members = $member_class->search({ $self->_foreign_key => $l_id }, { rows => 1 });
                unless ($members) {
                    $empty = $class->get(id => $l_id);
                    last;
                }
            }
        }
        $empty ||= $self->create($self->_member_key . 's' => []);
        
        $empty_lists{$class} = $empty;
        return $empty;
    }
    
    method _get_list (ClassName|Object $self: $orig, Maybe[Persistent] $id?, Maybe[ArrayRefOfPersistent] $members?) {
        $self->throw("You cannot supply both id and members") if $id && $members;
        
        if ($id) {
            return $self->$orig(id => $id);
        }
        elsif ($members) {
            if (@$members == 0) {
                return $self->_empty_list($orig);
            }
            
            my $member_class = $self->_member_class;
            
            my %lids;
            my $index = 0;
            foreach my $member (@$members) {
                ++$index;
                foreach my $lid ($member_class->get_column_values($self->_foreign_key, { $self->_member_key => $member->id })) {
                    $lids{$lid}++;
                }
            }
            
            my $expected_count = @$members;
            my $foreign_key    = $self->_foreign_key;
            foreach my $lid (sort { $a <=> $b } keys %lids) {
                next unless $lids{$lid} == $expected_count;
                next unless $member_class->search({ $foreign_key => $lid }) == $expected_count;
                return $self->$orig(id => $lid);
            }
            return $self->create($self->_member_key . 's' => $members);
        }
        else {
            $self->throw("List->get requires id or members");
        }
    }
    
    method _create_list (ClassName|Object $self: $orig, ArrayRefOfPersistent $members) {
        # create a new row, then use the new id to create new ListMember rows
        # for each supplied member
        my $list = $self->$orig();
        
        my $index = 0;
        my @lm_args;
        my $foreign_key = $self->_foreign_key;
        my $member_key  = $self->_member_key;
        foreach my $member (@$members) {
            push(@lm_args, { $foreign_key => $list->id, $member_key => $member->id });
        }
        my $member_class = $self->_member_class;
        $member_class->bulk_create_or_update(@lm_args);
        
        return $list;
    }
    
    method _instantiated_members {
        my @return;
        my $member_key = $self->_member_key;
        foreach my $member ($self->members) {
            push(@return, $member->$member_key);
        }
        return @return;
    }
}

1;
