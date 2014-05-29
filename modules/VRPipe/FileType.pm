
=head1 NAME

VRPipe::FileType - a factory for FileTypes

=head1 SYNOPSIS
    
    my $filetype = VRPipe::FileType->create($type, { file => $path });

=head1 DESCRIPTION

A FileType provides some essential infomation and methods for dealing with
files of that type without having to know about the format of those files
yourself.

For example, you can find out how many header lines are in a file using
C<num_header_lines()>, without knowing anything about the file, as long as its
type has been defined (and corresponds to one of the FileTypes).

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011,2014 Genome Research Limited.

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

package VRPipe::FileType;
use VRPipe::Base::AbstractFactory;

implementation_does qw/VRPipe::FileTypeRole/;

sub create {
    my $self = shift;
    my $return;
    eval { $return = $self->SUPER::create(@_); };
    if ($@) {
        # allow arbitrary filetypes by falling back to 'any' type
        if ($@ =~ qr/Invalid implementation class|perhaps you forgot to load/) {
            my $arbitrary = shift;
            die $@ if ($arbitrary eq 'any' || length($arbitrary) > 4);
            $return = $self->SUPER::create('any', @_);
            
            #*** type() is read-only and I can't figure out how to make the any
            # class rw via inheritance, so we just directly hack the hash
            $return->{type} = $arbitrary;
        }
        else {
            die $@;
        }
    }
    
    return $return;
}

1;
