
=head1 NAME

VRPipe::Steps::quantisnp_reformat_gs_export - a step

=head1 DESCRIPTION

Converts the Genome Studio csv files into a format that is suitable for running
with the quantiSNP scripts.

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

class VRPipe::Steps::quantisnp_reformat_gs_export with VRPipe::StepRole {
    method options_definition {
        return {
            reformat_script   => VRPipe::StepOption->create(description => 'full path to convert_fcr_to_quanti-per_sample.pl',                       optional => 0),
            manifest_file     => VRPipe::StepOption->create(description => 'full path to manifest file HumanCoreExome-12v1-0_A.csv',                 optional => 0),
            cnv_analysis_type => VRPipe::StepOption->create(description => 'type of cnv analysis, added to file metadata for downstream processing', optional => 1, default_value => 'quantisnp'),
        };
    }
    
    method inputs_definition {
        return {
            stepOne_file_input_GS_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'Genome Studio file containing genotyping data for reformatting',
                metadata    => { sample => 'sample name for cell line', storage_path => 'full path to iRODS file', analysis_uuid => 'analysis_uuid' },
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self              = shift;
            my $options           = $self->options;
            my $reformat_script   = $options->{reformat_script};
            my $manifest_file     = $options->{manifest_file};
            my $cnv_analysis_type = $options->{cnv_analysis_type};
            my $req               = $self->new_requirements(memory => 500, time => 1);
            foreach my $gs_file (@{ $self->inputs->{stepOne_file_input_GS_file} }) {
                my $gs_path  = $gs_file->path;
                my $basename = $gs_file->basename . '.reformat';
                my $new_meta = { cnv_analysis_type => $cnv_analysis_type };
                $gs_file->add_metadata($new_meta);
                my $reformatted_file = $self->output_file(output_key => 'stepOne_file_output_reformatted_file', basename => "$basename", type => 'txt', metadata => $gs_file->metadata);
                my $out_path         = $reformatted_file->path;
                my $cmd_line         = "$reformat_script $gs_path $manifest_file > $out_path";
                $self->dispatch([$cmd_line, $req]);
            }
        };
    }
    
    method outputs_definition {
        return {
            stepOne_file_output_reformatted_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Reformatted file',
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
