=head1 NAME

VRPipe::FileType::loc - loc filetype

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

When the B<VRPipe> local scheduler is used to schedule and run jobs, it produces
messages when the job completes. This filetype is for those report files.

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

class VRPipe::FileType::loc extends VRPipe::FileType::txt {
    method _build_read_backwards {
        return 1;
    }
}

1;