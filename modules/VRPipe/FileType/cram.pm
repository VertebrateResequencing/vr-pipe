
=head1 NAME

VRPipe::FileType::cram - cram filetype

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

The CRAM format is described here:
L<http://www.ebi.ac.uk/ena/about/cram_toolkit> It is used to hold aligned
sequence data in a highly compressed form.

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012,2015 Genome Research Limited.

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

class VRPipe::FileType::cram extends VRPipe::FileType::bam {
    method check_type {
        return ($self->hts_file_type =~ /^CRAM/) ? 1 : 0;
    }

}

1;
