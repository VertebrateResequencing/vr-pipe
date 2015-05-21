
=head1 NAME

VRPipe::FileProtocol - a factory for FileProtocols

=head1 SYNOPSIS
    
    my $protocol = 'irods';
    my $fileprotocol = VRPipe::FileProtocol->create($protocol, { file => $path });
    
    my $cat = $fileprotocol->cat;
    my $cmd = "$cat | myexe -";

=head1 DESCRIPTION

A FileProtocol provides some useful methods for dealing with files stored under
certain protocols without having to know how that protocol behaves yourself.

For example, you can build a command line to 'cat' a file and pipe it in to
something, regardless of if the file is on local disk or in irods etc, as long
as its protocol corresponds to one of the FileProtocols.

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

package VRPipe::FileProtocol;
use VRPipe::Base::AbstractFactory;

implementation_does qw/VRPipe::FileProtocolRole/;

1;
