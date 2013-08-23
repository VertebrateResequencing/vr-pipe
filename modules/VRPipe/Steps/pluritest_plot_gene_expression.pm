
=head1 NAME

VRPipe::Steps::pluritest_plot_gene_expression - a step

=head1 DESCRIPTION

Uses pluriTest R package (pluriTest_commandLine_vrpipe.r) to produce plots and
ancillary data from reformatted Genome Studio  gene expression data to
determine the pluripotency of different stem cell lines

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>.

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

class VRPipe::Steps::pluritest_plot_gene_expression with VRPipe::StepRole {
    method options_definition {
        return {
            pluritest_script => VRPipe::StepOption->create(description => 'path to modified pluritest R script'),
            pluritest_data   => VRPipe::StepOption->create(description => 'path to RData required by R script to plot graphs'),
            r_bin_path       => VRPipe::StepOption->create(description => 'path to specific version of R (2.15) required to run pluritest R script'),
        };
    }
    
    method inputs_definition {
        return {
            conv_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Files with reformatted Genom Studio gene expression data',
                max_files   => -1,
                min_files   => 0,
                metadata    => { merge_tag_id => 'tag id to enable sample to be identified in the multi-sample profile file', lanes => 'comma-separated list of lanes that the pluritest analysis is being performed on' },
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self             = shift;
            my $options          = $self->options;
            my $req              = $self->new_requirements(memory => 500, time => 1);
            my $pluritest_script = $options->{pluritest_script};
            my $pluritest_data   = $options->{pluritest_data};
            my $r_bin_path       = $options->{r_bin_path};
            
            foreach my $conv_file (@{ $self->inputs->{conv_files} }) {
                my $meta      = $conv_file->metadata;
                my $conv_path = $conv_file->path;
                $self->output_file(temporary  => 1,                 basename => 'pluritest_log.txt',      type => 'txt');
                $self->output_file(output_key => 'csv_file',        basename => 'pluritest.csv',          type => 'txt', metadata => $conv_file->metadata);
                $self->output_file(output_key => 'pluritest_plots', basename => 'pluritest_image01.png',  type => 'bin', metadata => $conv_file->metadata); #may need to add metadata
                $self->output_file(output_key => 'pluritest_plots', basename => 'pluritest_image02a.png', type => 'bin', metadata => $conv_file->metadata); #may need to add metadata
                $self->output_file(output_key => 'pluritest_plots', basename => 'pluritest_image02.png',  type => 'bin', metadata => $conv_file->metadata); #may need to add metadata
                $self->output_file(output_key => 'pluritest_plots', basename => 'pluritest_image03.png',  type => 'bin', metadata => $conv_file->metadata); #may need to add metadata
                $self->output_file(output_key => 'pluritest_plots', basename => 'pluritest_image03c.png', type => 'bin', metadata => $conv_file->metadata); #may need to add metadata
                my $cmd = qq[$r_bin_path --slave --args $conv_path $pluritest_data < $pluritest_script];
                $self->dispatch([$cmd, $req]);
            }
        };
    }
    
    method outputs_definition {
        return {
            pluritest_plots => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'png files produced by pluritest R script to assist in determination of pluripotency of stem cell lines',
                min_files   => 5,
                max_files   => -1,
                metadata    => { merge_tag_id => 'tag id to enable sample to be identified in the multi-sample profile file', lanes => 'comma-separated list of lanes that the pluritest analysis is being performed on' },
            ),
            csv_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'csv file with ancillary pluritest data per sample',
                max_files   => -1,
                metadata    => { merge_tag_id => 'tag id to enable sample to be identified in the multi-sample profile file', lanes => 'comma-separated list of lanes that the pluritest analysis is being performed on' },
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Uses R script (pluriTest_commandLine_vrpipe.r) to produce plots from reformatted Genome Studio gene expression data to determine pluripotency of stem cell lines";
    }
    
    method max_simultaneous {
        return 0;
    }
}

1;
