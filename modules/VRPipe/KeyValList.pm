
=head1 NAME

VRPipe::KeyValList - store an unordered list of VRPipe::KeyVal objects

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

You probably want to always use get() to safely create or get KeyValLists.

Direct use of create() will always create and return a new KeyValList, even if
one with the same set of keyvals already exists.

These lists are immuatble.

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

class VRPipe::KeyValList extends VRPipe::Persistent {
    __PACKAGE__->make_persistent(has_many => [members => 'VRPipe::KeyValListMember']);
    
    our $empty_kvl;
    
    around get (ClassName|Object $self: Persistent :$id?, ArrayRef[VRPipe::KeyVal] :$keyvals?, HashRef :$hash?) {
        $self->throw("You cannot supply both id and keyvals/hash") if ($id && ($keyvals || $hash));
        
        if ($id) {
            return $self->$orig(id => $id);
        }
        
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
        
        if ($keyvals) {
            if (@$keyvals == 0) {
                # we need to return the same keyvallist every time, but the
                # below code doesn't find those with no members, so we need to
                # special-case
                unless ($empty_kvl) {
                    #*** I have no idea how to do this in DBIx::Class...
                    # SELECT *
                    # FROM keyvallist
                    # LEFT OUTER JOIN keyvallistmember
                    #   ON (keyvallistmember.keyvallist = keyvallist.id)
                    #   WHERE keyvallistmember.keyvallist IS NULL
                    my $pager = VRPipe::KeyValList->get_column_values_paged('id', {});
                    PAGES: while (my $ids = $pager->next) {
                        foreach my $kvl_id (@$ids) {
                            my $members = VRPipe::KeyValListMember->search({ keyvallist => $kvl_id }, { rows => 1 });
                            unless ($members) {
                                $empty_kvl = VRPipe::KeyValList->get(id => $kvl_id);
                                last PAGES;
                            }
                        }
                    }
                }
                return $empty_kvl if $empty_kvl;
            }
            
            my %kvlids;
            my $index = 0;
            foreach my $keyval (@$keyvals) {
                ++$index;
                foreach my $kvlid (VRPipe::KeyValListMember->get_column_values('keyvallist', { keyval => $keyval->id })) {
                    $kvlids{$kvlid}++;
                }
            }
            
            my $expected_count = @$keyvals;
            foreach my $kvlid (sort { $a <=> $b } keys %kvlids) {
                next unless $kvlids{$kvlid} == $expected_count;
                next unless VRPipe::KeyValListMember->search({ keyvallist => $kvlid }) == $expected_count;
                return $self->$orig(id => $kvlid);
            }
            return $self->create(keyvals => $keyvals);
        }
        
        $self->throw("KeyValList->get requires id or keyvals or hash");
    }
    
    around create (ClassName|Object $self: ArrayRef[VRPipe::KeyVal] :$keyvals!) {
        # create a new row, then use the new id to create new
        # KeyValList rows for each supplied keyval
        my $list = $self->$orig();
        
        my $index = 0;
        my @kvm_args;
        foreach my $keyval (@$keyvals) {
            push(@kvm_args, { keyvallist => $list->id, keyval => $keyval->id });
        }
        VRPipe::KeyValListMember->bulk_create_or_update(@kvm_args);
        
        return $list;
    }
    
    method keyvals {
        my @return;
        foreach my $member ($self->members) {
            push(@return, $member->keyval);
        }
        return @return;
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
