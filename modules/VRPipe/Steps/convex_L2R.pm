
=head1 NAME

VRPipe::Steps::convex_L2R - a step

=head1 DESCRIPTION

Runs the SampleLogRatio R script from the CoNVex packages. Runs once per
pipeline, generating a log2 ratio file for each read depth file and a single
features file and correlation matrix file.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012-2014 Genome Research Limited.

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

class VRPipe::Steps::convex_L2R extends VRPipe::Steps::r_script {
    around options_definition {
        return {
            %{ $self->$orig },
            'regions_file'        => VRPipe::StepOption->create(description => 'regions file for which to generate read depths'),
            'convex_rscript_path' => VRPipe::StepOption->create(description => 'full path to CoNVex R scripts'),
            'includeChrX'         => VRPipe::StepOption->create(description => 'indicates whether to include Chr X in calulation', optional => 1, default_value => 0),
            'version'             => VRPipe::StepOption->create(description => 'SampleLogRatio script version', optional => 1, default_value => 3),
            'rpkm'                => VRPipe::StepOption->create(description => 'Use RPKM (aka FPKM) to correlate the samples', optional => 1, default_value => 0),
            'minSamples'          => VRPipe::StepOption->create(description => 'Minimum #samples expected in the subset when estimating median reference', optional => 1, default_value => 25),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => 'bam files',                                                metadata => { sample => 'sample name', sex => 'sample sex' }),
            rd_files  => VRPipe::StepIODefinition->create(type => 'rd',  max_files => -1, description => 'Set of convex Read Depths files (eg grouped by metadata)', metadata => { sample => 'sample name' }),
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $regions_file        = $options->{'regions_file'};
            my $convex_rscript_path = $options->{'convex_rscript_path'};
            my $includeChrX         = $options->{'includeChrX'};
            my $version             = $options->{'version'};
            my $rpkm                = $options->{'rpkm'};
            my $minSamples          = $options->{'minSamples'};
            
            my %samples;
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $meta   = $bam->metadata;
                my $sample = $meta->{sample};
                $samples{$sample}{bam} = $bam;
                $samples{$sample}{sex} = $meta->{sex};
            }
            foreach my $txt (@{ $self->inputs->{rd_files} }) {
                my $sample = $txt->metadata->{sample};
                $samples{$sample}{rd} = $txt if ($txt->type eq 'rd');
            }
            
            my $req = $self->new_requirements(memory => 2000, time => 1);
            
            # Create a tab-seperated SampleInfo file from the read depth for SampleLogRatioCall.R
            my $sample_info = $self->output_file(output_key => 'sample_info_file', basename => "sample_info.txt", type => 'txt');
            my $ofh = $sample_info->openw;
            
            foreach my $sample (keys %samples) {
                my $rd_file  = $samples{$sample}{rd};
                my $bam_path = $samples{$sample}{bam};
                my $sex      = $samples{$sample}{sex};
                
                my $basename = $rd_file->basename;
                $basename =~ s/\.rd$/.l2r/;
                
                my $rd_path = $rd_file->path;
                
                my $l2r_file = $self->output_file(output_key => 'l2r_files', basename => $basename, output_dir => $rd_file->dir, type => 'l2r', metadata => $rd_file->metadata);
                my $l2r_path = $l2r_file->path;
                
                print $ofh join("\t", $sample, $sex, $bam_path, $rd_path, $l2r_path), "\n";
            }
            $sample_info->close;
            my $sample_info_path = $sample_info->path;
            
            my $features_file = $self->output_file(output_key => 'features_file', basename => 'features.fts', type => 'fts');
            my $features_file_path = $features_file->path;
            
            my $corr_matrix_file = $self->output_file(output_key => 'corr_matrix_file', basename => 'corr_matrix.corr', type => 'corr');
            my $corr_matrix_file_path = $corr_matrix_file->path;
            
            my $sample_means_file = $self->output_file(output_key => 'sample_means_file', basename => 'sample_means.savg', type => 'savg');
            my $sample_means_file_path = $sample_means_file->path;
            
            my $cmd = $self->rscript_cmd_prefix . " $convex_rscript_path/SampleLogRatioCall.R $sample_info_path,$regions_file,$features_file_path,$includeChrX,$corr_matrix_file_path,$version,$rpkm,$minSamples,$sample_means_file_path";
            $self->dispatch([$cmd, $req]);
        };
    }
    
    method outputs_definition {
        return {
            sample_info_file  => VRPipe::StepIODefinition->create(type => 'txt',  max_files => 1,  description => 'sample info file listing sample, sex, bam_path, read_depth_path, l2r_path'),
            features_file     => VRPipe::StepIODefinition->create(type => 'fts',  max_files => 1,  description => 'a single convex features file for each set of L2R files'),
            corr_matrix_file  => VRPipe::StepIODefinition->create(type => 'corr', max_files => 1,  description => 'a single convex correlation matrix file for each set of L2R files'),
            sample_means_file => VRPipe::StepIODefinition->create(type => 'savg', max_files => 1,  description => 'a single convex sample means file'),
            l2r_files         => VRPipe::StepIODefinition->create(type => 'l2r',  max_files => -1, description => 'a log 2 ratio file for each input read depths file', metadata => { sample => 'sample name' }),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs CoNVex SampleLogRatioCall, generating log 2 ratio files from a fofn of read depths files and a features file and correlation matrix file and a sample means file";
    }
    
    method max_simultaneous {
        return 1;            # should only run once
    }
}

1;
