
=head1 NAME

VRPipe::Steps::reformat_cnv_output_to_bed - a step

=head1 DESCRIPTION

reformat quantisnp or penncnv cnv output files to bed format.

=head1 AUTHOR

Phil Carter <pc12@sanger.ac.uk>, John Maslen <jm23@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::reformat_cnv_output_to_bed with VRPipe::StepRole {
    method options_definition {
        return {
            control_metadata_key => VRPipe::StepOption->create(description => 'metadata key that determines the control sample (1=control, 0=query sample)', optional => 1, default_value => 'sample_control'),
            cnv_analysis_type    => VRPipe::StepOption->create(description => 'Determines the type of analysis, penncnv or quantisnp',                       optional => 1, default_value => 'penncnv'),
        };
    }
    
    method inputs_definition {
        return {
            cnv_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'CNV output files from quantisnp or penncnv pipelines',
                metadata    => { sample => 'sample name', cnv_analysis_type => 'quantisnp or penncnv' },
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self              = shift;
            my $options           = $self->options;
            my $control_key       = $options->{control_metadata_key};
            my $cnv_analysis_type = $options->{cnv_analysis_type};
            my $req               = $self->new_requirements(memory => 500, time => 1);
            foreach my $cnv_file (@{ $self->inputs->{cnv_files} }) {
                my $cnv_path = $cnv_file->path;
                my $meta     = $cnv_file->metadata;
                my $sample   = $meta->{sample};
                my $basename = $meta->{$control_key} ? $sample . '_' . $cnv_analysis_type . '_CONTROL.bed' : $sample . '_' . $cnv_analysis_type . '.bed';
                my $bed_file = $self->output_file(output_key => 'bed_files', basename => "$basename", type => 'txt', metadata => $cnv_file->metadata);
                my $bed_path = $bed_file->path;
                my $cmd      = "use VRPipe::Steps::reformat_cnv_output_to_bed; VRPipe::Steps::reformat_cnv_output_to_bed->bed_conversion_output(q[$cnv_path], q[$bed_path], q[$cnv_analysis_type]);";
                $self->dispatch_vrpipecode($cmd, $req, { output_files => [$bed_file] });
            }
        };
    }
    
    method outputs_definition {
        return {
            bed_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Converted bed format file for CNV comparison',
                max_files   => -1,
                min_files   => 0,
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Reformats CNV output file to bed format";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method bed_conversion_output (ClassName|Object $self: Str|File $in_cnv!, Str|File $out_bed!, Str $cnv_type!) {
        my $bed_file = VRPipe::File->get(path => $out_bed);
        
        my $bed_fh = $bed_file->openw;
        
        my $expected_lines = 0;
        
        my $fh;
        open($fh, $in_cnv) || $self->throw("Couldn't open '$in_cnv': $!");
        
        while (<$fh>) {
            chomp;
            my @cnv_chr_arr = (split /\s+/);
            if ($cnv_type eq 'penncnv') {
                @cnv_chr_arr = split(/^chr(.*)\:([0-9]+)\-([0-9]+)/, $cnv_chr_arr[0]);
            }
            if ($cnv_chr_arr[3] =~ /^[0-9]+$/) {
                print $bed_fh join("\t", @cnv_chr_arr[1 .. 3]), "\n";
                $expected_lines++;
            }
        }
        $bed_file->close;
        
        my $actual_lines = $bed_file->lines;
        if ($actual_lines == $expected_lines) {
            return 1;
        }
        else {
            $bed_file->unlink;
            $self->warn("Wrote $expected_lines lines to " . $bed_file->path . ", but only read back $actual_lines! Deleted the output.");
            return 0;
        }
    }
}

1;
