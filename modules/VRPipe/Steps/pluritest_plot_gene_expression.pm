
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

class VRPipe::Steps::pluritest_plot_gene_expression extends VRPipe::Steps::r {
    around options_definition {
        return {
            %{ $self->$orig },
            pluritest_script => VRPipe::StepOption->create(description => 'path to modified pluritest R script'),
            pluritest_data   => VRPipe::StepOption->create(description => 'path to RData required by R script to plot graphs')
        };
    }
    
    method inputs_definition {
        return {
            conv_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'File with reformatted Genome Studio gene expression data',
                min_files   => 0,
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            my $req              = $self->new_requirements(memory => 500, time => 1);
            my $pluritest_script = $options->{pluritest_script};
            my $pluritest_data   = $options->{pluritest_data};
            
            my ($conv_file) = @{ $self->inputs->{conv_files} };
            my $meta        = $conv_file->metadata;
            my $conv_path   = $conv_file->path;
            $self->output_file(temporary  => 1,          basename => 'pluritest_log.txt', type => 'txt', metadata => $meta);
            $self->output_file(output_key => 'csv_file', basename => 'pluritest.csv',     type => 'txt', metadata => $meta);
            foreach my $basename (qw(pluritest_image01.png pluritest_image02a.png pluritest_image02.png pluritest_image03.png pluritest_image03c.png)) {
                $self->output_file(output_key => 'pluritest_plots', basename => $basename, type => 'bin', metadata => $meta);
            }
            
            my $cmd = $self->r_cmd_prefix . qq[ --slave --args $conv_path $pluritest_data < $pluritest_script];
            $self->dispatch([$cmd, $req]);
        };
    }
    
    method outputs_definition {
        return {
            pluritest_plots => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'png files produced by pluritest R script to assist in determination of pluripotency of stem cell lines',
                min_files   => 5,
                max_files   => -1,
            ),
            csv_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'csv file with ancillary pluritest data per sample',
                max_files   => -1,
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
