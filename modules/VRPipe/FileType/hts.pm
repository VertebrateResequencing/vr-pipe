
=head1 NAME

VRPipe::FileType::hts - hts filetype base class

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

class VRPipe::FileType::hts with VRPipe::FileTypeRole {
    my $htsfile_exe = file($ENV{HTSLIB}, 'bin', 'htsfile');
    
    method hts_file_type {
        my $path = $self->file;
        return `$htsfile_exe $path | cut -f2`;
    }
    
    method num_records {
        my $path = $self->file;
        open(my $wc, "$htsfile_exe -cH $path | wc -l |") || $self->throw("$htsfile_exe -cH $path | wc -l did not work");
        my ($lines) = split(" ", <$wc>);
        CORE::close($wc);
        return $lines;
    }
    
    method num_header_lines {
        return scalar(@{ $self->header_lines });
    }
    
    method header_lines {
        my $path         = $self->file;
        my $headers      = `$htsfile_exe -ch $path`;
        my @header_lines = split(/\n/, $headers);
        return \@header_lines;
    }
}

1;
