
=head1 NAME

VRPipe::FileType::fq - fastq filetype

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

Fastq files (often abbrievated to 'fq') hold sequence and quality information
produced by genome sequencers. The format is described here:
L<http://maq.sourceforge.net/fastq.shtml>

*** more documentation to come

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

class VRPipe::FileType::fq extends VRPipe::FileType::txt {
    #*** needs more implementation!
    
    method check_type {
        my $file = $self->file;
        my $type = $self->type;
        $file =~ s/\.gz$// unless $type eq 'gz';
        if ($file =~ /\.(?:$type|fastq)$/) { #*** this sucks as a test...
            return 1;
        }
        return 0;
    }
    
    method num_records {
        my $path = $self->file;
        my $records;
        if ($path =~ /\.gz$/) {
            $records = `gunzip -c $path | fastqcheck | head -1`;
        }
        else {
            $records = `fastqcheck $path | head -1`;
        }
        ($records) = $records =~ m/^(\d+) sequences/;
        $records |= 0;
        return $records;
    }
}

1;
