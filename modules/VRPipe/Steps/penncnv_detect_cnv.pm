
=head1 NAME

VRPipe::Steps::penncnv_detect_cnv - a step

=head1 DESCRIPTION

Detects raw CNVs using the Genome Studio genotyping files

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

class VRPipe::Steps::penncnv_detect_cnv with VRPipe::StepRole {
    method options_definition {
        return {
            detect_cnv_script => VRPipe::StepOption->create(description => 'full path to detect_cnv.pl',                                             optional => 0),
            detect_cnv_hmm    => VRPipe::StepOption->create(description => 'full path to custom.hmm',                                                optional => 0),
            detect_cnv_pfb    => VRPipe::StepOption->create(description => 'full path to PFB',                                                       optional => 0),
            cnv_analysis_type => VRPipe::StepOption->create(description => 'type of cnv analysis, added to file metadata for downstream processing', optional => 1, default_value => 'penncnv'),
        };
    }
    
    method inputs_definition {
        return {
            stepOne_file_input_GS_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'Genome Studio file containing genotyping data for calling CNVs',
                metadata    => { sample => 'sample name for cell line', storage_path => 'full path to iRODS file', analysis_uuid => 'analysis_uuid' },
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self               = shift;
            my $options            = $self->options;
            my $detect_cnv_script  = $options->{detect_cnv_script};
            my $detect_cnv_hmm     = $options->{detect_cnv_hmm};
            my $detect_cnv_pfb     = $options->{detect_cnv_pfb};
            my $cnv_analysis_type  = $options->{cnv_analysis_type};
            my $detect_cnv_options = "-test -hmm $detect_cnv_hmm -pfb $detect_cnv_pfb --confidence --region 1-22 --minsnp 1";
            my $req                = $self->new_requirements(memory => 500, time => 1);
            foreach my $gs_file (@{ $self->inputs->{stepOne_file_input_GS_file} }) {
                my $gs_path  = $gs_file->path;
                my $basename = $gs_file->basename . '.rawcnv';
                my $new_meta = { cnv_analysis_type => $cnv_analysis_type };
                $gs_file->add_metadata($new_meta);
                my $raw_cnv_file = $self->output_file(output_key => 'stepOne_file_output_raw_cnv_file', basename => "$basename", type => 'txt', metadata => $gs_file->metadata);
                my $out_path     = $raw_cnv_file->path;
                my $logFile      = $out_path . '.LOG';
                my $cmd_line     = "perl $detect_cnv_script $detect_cnv_options -log $logFile -out $out_path $gs_path";
                $self->dispatch([$cmd_line, $req]);
            }
        };
    }
    
    method outputs_definition {
        return {
            stepOne_file_output_raw_ncv_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Raw CNV data',
                max_files   => -1,                                                                                                                    # -1 = As many as you like
                min_files   => 0,
                metadata    => { sample => 'sample name for cell line', storage_path => 'full path to iRODS file', analysis_uuid => 'analysis_uuid' },
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Detects raw CNVs using the Genome Studio genotyping files";
    }
    
    method max_simultaneous {
        return 0;
    }

}

1;
