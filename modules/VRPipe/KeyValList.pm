
=head1 NAME

VRPipe::KeyValList - store an unordered list of VRPipe::KeyVal objects

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

These lists are immutable and insensitive to order.

Both get() and create() will return an existing list if one had previously been
created using the supplied list of members. Both will also create and return a
new list if one had not been previously created.

NB: This looks like a normal PersistentList class, but its member instances
directly hold the desired information, instead of referencing a keyval class.

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
    use VRPipe::Persistent::InMemory;
    
    sub _member_class { 'VRPipe::KeyValListMember' }
    sub _member_key   { 'keyval' }
    sub _foreign_key  { 'keyvallist' }
    
    has 'lookup' => (
        is     => 'rw',
        isa    => Varchar [64],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    __PACKAGE__->make_persistent(has_many => [members => 'VRPipe::KeyValListMember']);
    
    sub _members_to_string {
        my ($self, $members) = @_;
        return 'undef' unless @$members;
        my @keyvals = sort map { $_->[0] . (ref($_->[1]) eq 'ARRAY' ? join(',', @{ $_->[1] }) : $_->[1]) } @$members;
        return join('', @keyvals);
    }
    
    around get (ClassName|Object $self: Persistent :$id?, ArrayRef :$keyvals?, HashRef :$hash?) {
        if ($hash) {
            $keyvals ||= [];
            while (my ($key, $val) = each %{$hash}) {
                push(@$keyvals, [$key, $val]);
            }
        }
        
        return $self->_get_list($orig, $id, $keyvals);
    }
    
    around create (ClassName|Object $self: ArrayRef :$keyvals!) {
        # because we don't have a real keyval class but instead store data on
        # the Member class directly for speed reasons, we have to implement this
        # ourselves
        my $lookup = $self->_string_to_lookup($self->_members_to_string($keyvals));
        
        my $im       = VRPipe::Persistent::InMemory->new();
        my $lock_key = 'PersistentList.' . $lookup;
        $im->block_until_locked($lock_key);
        
        my $list = $self->$orig(lookup => $lookup);
        
        my @lm_args;
        my $lid = $list->id;
        foreach my $ref (@$keyvals) {
            my ($key, $val) = @$ref;
            my @vals = ref($val) eq 'ARRAY' ? @$val : ($val);
            foreach my $value (@vals) {
                push(@lm_args, { keyvallist => $lid, keyval_key => $key, val => $value });
            }
        }
        VRPipe::KeyValListMember->bulk_create_or_update(@lm_args);
        
        $im->unlock($lock_key);
        return $list;
    }
    
    method keyvals {
        # (_instantiated_members assumes KeyValListMember has a column that
        # points to another table, but we don't have that for speed reasons)
        my $ref = VRPipe::KeyValListMember->get_column_values(['keyval_key', 'val'], { keyvallist => $self->id });
        return @$ref;
    }
    
    # speed critical, so sub instead of method, and we reimplement keyvals to
    # avoid the call stack. Note, if the list had the same key multiple times,
    # our returned hash only contains the last value for that key
    sub as_hashref {
        my $self   = shift;
        my $return = {};
        my $ref    = VRPipe::KeyValListMember->get_column_values(['keyval_key', 'val'], { keyvallist => $self->id });
        foreach my $ref (@$ref) {
            my ($keyval_key, $val) = @$ref;
            if (exists $return->{$keyval_key}) {
                my $previous = $return->{$keyval_key};
                if (ref($previous) eq 'ARRAY') {
                    push(@{ $return->{$keyval_key} }, $val);
                }
                else {
                    $return->{$keyval_key} = [$previous, $val];
                }
            }
            else {
                $return->{$keyval_key} = $val;
            }
        }
        return $return;
    }
    
    method get_value (Str $key) {
        my @vals = VRPipe::KeyValListMember->get_column_values('val', { keyvallist => $self->id, keyval_key => $key });
        if (wantarray) {
            return @vals;
        }
        else {
            return pop(@vals);
        }
    }
    
    # we don't have an add_value or similar, because the list must be immutable;
    # they are used to store metadata on File rows, and 2 different files could
    # have the same metadata and therefore store the same keyvallist foreign
    # key id. So we can't just alter the keyvallist when the metadata on one
    # of the files changes, because that would affect the other as well
}

1;
