
=head1 NAME

VRPipe::FileProtocol::irods - handler for files stored in irods

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

Handles irods files with knowledge of how to get their file content streamed
out for piping purposes.

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

class VRPipe::FileProtocol::irods with VRPipe::FileProtocolRole {
    method cat {
        my $file = $self->file;
        my $cat  = "iget $file -";
        if ($file =~ /\.gz$/) {
            $cat .= ' | gzip -dc';
        }
        return $cat;
    }
}

1;
