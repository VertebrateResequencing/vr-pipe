
=head1 NAME

VRPipe::Steps::convex_L2R - a step

=head1 DESCRIPTION

Runs the SampleLogRatio R script from the CoNVex packages. Runs once per
pipeline, generating a log2 ratio file for each read depth file and a single
features file.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Steps::convex_L2R extends VRPipe::Steps::r_script {
    around options_definition {
        return { %{ $self->$orig },
                 'regions_file' => VRPipe::StepOption->create(description => 'regions file for which to generate read depths'),
                 'convex_rscript_path' => VRPipe::StepOption->create(description => 'full path to CoNVex R scripts'),
                 'includeChrX'  => VRPipe::StepOption->create(description => 'indicates whether to include Chr X in calulation', optional => 1, default_value => 1), };
    }
    
    method inputs_definition {
        return { rd_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'Set of convex Read Depths files (eg grouped by metadata)'), };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options      = $self->options;
            $self->handle_standard_options($options);

            my $regions_file = $options->{'regions_file'};
            my $convex_rscript_path = $options->{'convex_rscript_path'};
            my $includeChrX  = $options->{'includeChrX'};
            
            my $req = $self->new_requirements(memory => 2000, time => 1);
            
            # Create a tab-seperated SampleInfo file from the read depth for SampleLogRatioCall.R
            my $sample_info = $self->output_file(basename => "SampleInfo.txt", type => 'txt', temporary => 1);
            my $ofh = $sample_info->openw;
            
            my @l2r_files;
            
            foreach my $rd_file (@{ $self->inputs->{rd_files} }) {
                my $basename = $rd_file->basename;
                $basename =~ s/\.rd\.txt$/.l2r.txt/;
                
                my $rd_dir  = $rd_file->dir;
                my $rd_path = $rd_file->path;
                
                my $bam_path   = $rd_file->metadata->{source_bam};
                my $bam_file   = VRPipe::File->get(path => $bam_path);
                my $bam_sample = $bam_file->metadata->{sample};
                my $bam_sex    = $bam_file->metadata->{sex};
                
                my $l2r_file = $self->output_file(output_key => 'l2r_files', basename => $basename, output_dir => $rd_dir, type => 'txt');
                push(@l2r_files, $l2r_file);
                my $l2r_path = $l2r_file->path;
                
                print $ofh join("\t", $bam_sample, $bam_sex, $bam_path, $rd_path, $l2r_path), "\n";
            }
            $sample_info->close;
            my $sample_info_path = $sample_info->path;
            
            my $features_file = $self->output_file(output_key => 'features_file', basename => 'features.txt', type => 'txt');
            my $features_file_path = $features_file->path;
            push(@l2r_files, $features_file);
            
            my $cmd =  $self->rscript_cmd_prefix . " $convex_rscript_path/SampleLogRatioCall.R $sample_info_path,$regions_file,$features_file_path,$includeChrX";
            
            $self->dispatch([$cmd, $req, { output_files => \@l2r_files }]);
        
        };
    }
    
    method outputs_definition {
        return { features_file => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1,  description => 'a single convex features file for each set of L2R files'),
                 l2r_files     => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'a log 2 ratio file for each input read depths file'), };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs CoNVex SampleLogRatioCall, generating log 2 ratio files from a fofn of read depths files and a single features file";
    }
    
    method max_simultaneous {
        return 1;            # should only run once
    }
}

1;
