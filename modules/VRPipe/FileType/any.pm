
=head1 NAME

VRPipe::FileType::any - any filetype

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A generic filetype for any file at all, for use when the filetype of a file
isn't known. If anything tries to work with files of this type, the effect will
most likely be similar to treating the file as a text file.

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

use VRPipe::Base;

class VRPipe::FileType::any with VRPipe::FileTypeRole {
    method check_type {
        my $type = $self->type;
        return 1 if $type eq 'any';
        
        #*** can't get inheritance from FileTypeRole to work, so this is
        # duplicated code :(
        my $file = $self->file;
        $file =~ s/\.gz$// unless $type eq 'gz';
        if ($file =~ /\.$type$/) {
            return 1;
        }
        return 0;
    }
}

1;
