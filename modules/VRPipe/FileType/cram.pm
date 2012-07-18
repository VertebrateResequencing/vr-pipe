
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

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::FileType::cram extends VRPipe::FileType::bin {
    our $correct_magic = [qw(103 122 101 115 000 000 000 103 037 213 010 000 000 000 000 000)];
    
    around check_type {
        $self->$orig || return 0;
        my $path = $self->file;
        return $path =~ /\.cram$/ ? 1 : 0;
    }
    
    method num_header_lines (ClassName|Object $self: Str|File :$reference_fasta!) {
        my $path         = $self->file;
        my $headers      = `java -Dreference=$reference_fasta -cp $ENV{CRAMTOOLS}/cramtools.jar net.sf.picard.sam.ViewSam INPUT=$path | samtools view -SH -`;
        my @header_lines = split(/\n/, $headers);
        return scalar(@header_lines);
    }
    
    method num_records (ClassName|Object $self: Str|File :$reference_fasta!) {
        my $path    = $self->file;
        my $records = `java -Dreference=$reference_fasta -cp $ENV{CRAMTOOLS}/cramtools.jar net.sf.picard.sam.ViewSam INPUT=$path | samtools view -Sc -`;
        ($records) = $records =~ /^(\d+)/m;
        $records ||= 0;
        return $records;
    }
}

1;
