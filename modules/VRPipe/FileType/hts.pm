
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
}

1;
