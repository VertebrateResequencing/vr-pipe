
=head1 NAME

VRPipe::FileType::vcf - VCF filetype

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

The VCF format is described here: L<http://vcftools.sourceforge.net/specs.html>
It is used to describe variants in biological sequences.

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

class VRPipe::FileType::vcf extends VRPipe::FileType::hts {
    method check_type {
        return ($self->hts_file_type =~ /^VCF/) ? 1 : 0;
    }
    
    method samples {
        my $header_lines = $self->header_lines;
        $self->throw("No header lines found in, " . $self->file) unless (@$header_lines);
        my @line = split(/\t/, ${$header_lines}[-1]);
        my @samples = @line[9 .. $#line];
        @samples || $self->throw("No samples found in " . $self->file);
        return \@samples;
    }
    
    method num_samples {
        return scalar(@{ $self->samples });
    }

}

1;
