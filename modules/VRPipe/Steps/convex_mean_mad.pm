
=head1 NAME

VRPipe::Steps::convex_mean_mad - a step

=head1 DESCRIPTION

XXX

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::convex_mean_mad extends VRPipe::Steps::r_script {
    around options_definition {
        return {
            %{ $self->$orig },
            'regions_file'        => VRPipe::StepOption->create(description => 'regions file for which to generate read depths'),
            'features_file'       => VRPipe::StepOption->create(description => 'features file from L2R calculation step'),
            'corr_matrix_file'    => VRPipe::StepOption->create(description => 'correlation matrix file from L2R calculation step'),
            'convex_rscript_path' => VRPipe::StepOption->create(description => 'full path to CoNVex R scripts'),
            'convex_classpath'    => VRPipe::StepOption->create(description => 'convex classpath to all convex package jars'),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => 'bam files',                            metadata => { sample => 'sample name', sex => 'sample sex' }),
            rd_files  => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'read depth (rd) txt files per-sample', metadata => { sample => 'sample name' }),
            l2r_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'l2r txt files per-sample',             metadata => { sample => 'sample name' }),
            gam_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'gam txt files per-sample',             metadata => { sample => 'sample name' }),
            cnv_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'convex CNV txt files per-sample',      metadata => { sample => 'sample name' }),
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $regions_file        = $options->{'regions_file'};
            my $features_file       = $options->{'features_file'};
            my $corr_matrix_file    = $options->{'corr_matrix_file'};
            my $convex_rscript_path = $options->{'convex_rscript_path'};
            my $convex_classpath    = $options->{'convex_classpath'};
            
            my %samples;
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $meta   = $bam->metadata;
                my $sample = $meta->{sample};
                $samples{$sample}{bam} = $bam;
                $samples{$sample}{sex} = $meta->{sex};
            }
            foreach my $txt (@{ $self->inputs->{rd_files} }, @{ $self->inputs->{l2r_files} }, @{ $self->inputs->{gam_files} }, @{ $self->inputs->{cnv_files} }) {
                my $sample = $txt->metadata->{sample};
                $samples{$sample}{rd}  = $txt->path if ($txt->path =~ /rd.txt$/);
                $samples{$sample}{l2r} = $txt->path if ($txt->path =~ /l2r.txt$/);
                $samples{$sample}{gam} = $txt->path if ($txt->path =~ /gam.txt$/);
                if ($txt->path =~ /cnvs.txt$/) {
                    $samples{$sample}{cnvs}   = $txt->path;
                    $samples{$sample}{outdir} = $txt->dir;
                    my $base = $txt->basename;
                    $base =~ s/txt$/mm.txt/;
                    $samples{$sample}{basename} = $base;
                }
            }
            
            # Create a tab-seperated SampleInfo file from the read depth for SampleLogRatioCall.R
            my $sample_info = $self->output_file(output_key => 'sample_info_file', basename => "SampleInfo.txt", type => 'txt');
            my $ofh = $sample_info->openw;
            
            foreach my $sample (keys %samples) {
                my $bam       = $samples{$sample}{bam};
                my $bam_path  = $bam->path;
                my $sex       = $samples{$sample}{sex};
                my $rd_path   = $samples{$sample}{rd};
                my $l2r_path  = $samples{$sample}{l2r};
                my $gam_path  = $samples{$sample}{gam};
                my $cnvs_path = $samples{$sample}{cnvs};
                
                my $basename = $samples{$sample}{basename};
                my $outdir   = $samples{$sample}{outdir};
                
                my $mean_mad_file = $self->output_file(output_key => 'cnv_files_with_mean_mad', basename => $basename, output_dir => $outdir, type => 'txt', metadata => $bam->metadata);
                my $mean_mad_path = $mean_mad_file->path;
                
                $samples{$sample}{mean_mad} = $mean_mad_file;
                
                print $ofh join("\t", $sample, $sex, $bam_path, $rd_path, $l2r_path, $gam_path, $cnvs_path, $mean_mad_path), "\n";
            }
            $sample_info->close;
            my $sample_info_path = $sample_info->path;
            
            my $req = $self->new_requirements(memory => 2000, time => 1);
            foreach my $sample (keys %samples) {
                my $cnv_path      = $samples{$sample}{cnvs};
                my $mean_mad_file = $samples{$sample}{mean_mad};
                my $mean_mad_path = $mean_mad_file->path;
                my $cmd           = "export CLASSPATH $convex_classpath; " . $self->rscript_cmd_prefix . " $convex_rscript_path/MeansMadsCall.R $cnv_path,$sample,$corr_matrix_file,$sample_info_path,$regions_file,$features_file,.,$mean_mad_path";
                $self->dispatch([$cmd, $req, { output_files => [$mean_mad_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            sample_info_file        => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1,  description => 'text file listing sample, sex, bam_path, rd_path, l2r_path, gam_path, cnvs_path and mean_mad_path for all samples'),
            cnv_files_with_mean_mad => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'convex cnv calls with mean and mad statistics', metadata => { sample => 'sample name' }),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs CoNVex SampleLogRatioCall, generating log 2 ratio files from a fofn of read depths files and a features file and correlation matrix file";
    }
    
    method max_simultaneous {
        return 1;            # should only run once
    }
}

1;
