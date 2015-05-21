
=head1 NAME

VRPipe::FileType::bam - bam filetype

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

bam files are the compressed binary forms of sam files:
L<http://samtools.sourceforge.net/> They hold aligned biological sequence data.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::FileType::bam extends VRPipe::FileType::hts {
    my $samtools_exe = file($ENV{SAMTOOLS}, 'samtools');
    
    method check_type {
        return ($self->hts_file_type =~ /^BAM/) ? 1 : 0;
    }
    
    method num_records {
        my $path    = $self->file;
        my $sam     = -T $path ? 'S' : '';
        my $records = `$samtools_exe view -${sam}c -F 0x900 $path`;
        ($records) = $records =~ /^(\d+)/m;
        $records ||= 0;
        return $records;
    }
}

1;
