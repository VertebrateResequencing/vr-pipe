
=head1 NAME

VRPipe::PersistentLocklessCreate - Persistent with a lockless create method

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This is VRPipe:Persistent but with the create() method overridden so that
instead of doing a select FOR UPDATE, which can block up the database with
queries and take a long time to finish the select, it does a get first, and
only does the create on failure.

This is useful for classes that expect to create new rows rarely, and spend
most of their time getting existing rows via the create() call.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::PersistentLocklessCreate extends VRPipe::Persistent {
    sub create {
        my $self = shift;
        my $return = $self->_get(-1, @_);
        unless ($return) {
            # we really really don't want to resort to 'FOR UPDATE' queries in
            # MySQL, so we'll block and lock using Redis, which will be much
            # faster
            my $im        = VRPipe::Persistent::InMemory->new();
            my %args      = @_;
            my $fa        = $self->_find_args(\%args);
            my %find_args = %{ $fa->[0] || {} };
            my $key       = join(',', map { $_ . '=>' . $find_args{$_} } sort keys %find_args);
            if ($key) {
                $im->block_until_locked($key);
                $return = $self->_get(-1, @_);
            }
            unless ($return) {
                $return = $self->_get(1, @_);
            }
            $im->unlock($key) if $key;
        }
        return $return;
    }
}

1;
