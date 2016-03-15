
=head1 NAME

VRPipe::FileType::hts - hts filetype base class

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

Filetype detection is done using Inline C code that uses the htslib library,
which you need to have compiled and the main directory (containing include and
lib subdirectories) pointed to by the HTSLIB environment variable.

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>. Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015-2016 Genome Research Limited.

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
    
    use Inline C => Config => INC => "-I$ENV{HTSLIB}/include" => LIBS => "-L$ENV{HTSLIB}/lib -lhts -lz" => GLOBAL_LOAD => 1;
    
    method hts_file_type {
        my $format = $self->_c_hts_file_type($self->file->stringify);
        unless ($format) {
            $self->throw("Failed to cope with " . $self->file->stringify);
        }
        return $format;
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
    
    method check_records_vs_input (VRPipe::File $in_file, Str $cmd_line?) {
        my $out_file = VRPipe::File->get(path => $self->file);
        $out_file->update_stats_from_disc(retries => 3);
        
        my $reads         = $in_file->metadata->{reads};
        my $expected_recs = $reads || $in_file->num_records;
        my $actual_recs   = $out_file->num_records;
        $out_file->add_metadata({ reads => $actual_recs }) if $out_file->meta_value('reads');
        
        if ($actual_recs == $expected_recs) {
            return 1;
        }
        elsif ($reads && $reads > $actual_recs) {
            # metadata might be wrong or based on raw reads, not 0x900 reads, so
            # double-check
            $expected_recs = $in_file->num_records;
            if ($actual_recs == $expected_recs) {
                $in_file->add_metadata({ reads => $expected_recs });
                return 1;
            }
        }
        # else
        if ($cmd_line) {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_recs records were generated in the output file, yet there were $expected_recs records in the input file");
        }
        else {
            return 0;
        }
    }
    
    use Inline C => <<'END_C';

#include "htslib/hfile.h"
#include "htslib/hts.h"
#include <errno.h>
#include <stdio.h>

// code lifted from htslib htsfile.c authored by John Marshall

char* _c_hts_file_type(SV* self, char* path) {
    htsFormat fmt;
    hFILE *fp = hopen(path, "r");
    if (fp == NULL) {
        fprintf(stderr, "htsfile: can't open \"%s\": %s\n", path, strerror(errno));
        return 0;
    }
    
    if (hts_detect_format(fp, &fmt) < 0) {
        hclose_abruptly(fp);
        fprintf(stderr, "htsfile: detecting \"%s\" format failed: %s\n", path, strerror(errno));
        return 0;
    }
    
    char *description = hts_format_description(&fmt);
    return description;
}

END_C

}

1;
