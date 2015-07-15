
=head1 NAME

VRPipe::FileType::idx - generic filetype for hts index files: bai,crai,csi,tbi

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::FileType::idx extends VRPipe::FileType::hts {
    method check_type {
        my $type = $self->hts_file_type;
        if ($type =~ /^BAI/ || $type =~ /^CRAI/ || $type =~ /^CSI/ || $type =~ /^TBI/) {
            return 1;
        }
        elsif ($type =~ /^SAM/) {
            # this is a hack to get around this issue:
            # https://github.com/samtools/htslib/issues/200
            # which means CRAI files are detected as SAM files.
            # a fix for this is expected in the next htslib
            # release. in the meantime, we will be returning
            # true here for all gzip compressed text files,
            # which is not ideal
            my $file = $self->file;
            return -B $file ? 1 : 0;
        }
        else {
            return 0;
        }
    }

}

1;
