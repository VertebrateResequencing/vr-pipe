
=head1 NAME

VRPipe::Steps::pluritest_plot_gene_expression - a step

=head1 DESCRIPTION

Uses pluriTest R package (pluriTest_commandLine_vrpipe.r) to produce plots and
ancillary data from reformatted Genome Studio  gene expression data to
determine the pluripotency of different stem cell lines

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>. Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013-2015 Genome Research Limited.

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
            my $req              = $self->new_requirements(memory => 1000, time => 1);
            my $pluritest_script = $options->{pluritest_script};
            my $pluritest_data   = $options->{pluritest_data};
            
            my ($conv_file) = @{ $self->inputs->{conv_files} };
            my $meta        = $conv_file->metadata;
            my $conv_path   = $conv_file->path->stringify;
            $self->output_file(temporary => 1, basename => 'pluritest_log.txt', type => 'txt', metadata => $meta);
            my $csv_path = $self->output_file(output_key => 'csv_file', basename => 'pluritest.csv', type => 'txt', metadata => $meta)->path->stringify;
            $self->relate_input_to_output($conv_path, 'pluritest_summary', $csv_path);
            
            my $vrtrack = VRPipe::Schema->create("VRTrack");
            my $donor_node = $vrtrack->get('Donor', { id => $meta->{sample_cohort} });
            $self->throw("No donor node with id $meta->{sample_cohort}") unless $donor_node;
            
            # first remove any old pluritest plots that are attached to the donor
            foreach my $plot ($donor_node->related(outgoing => { type => 'pluritest_plot' })) {
                $donor_node->divorce_from($plot);
            }
            
            my %num_to_type = ('01' => 'intensity', '02a' => 'clustering', '02' => 'pluripotency', '03' => 'pluripotency_vs_novelty', '03c' => 'novelty');
            foreach my $num (qw(01 02a 02 03 03c)) {
                my $basename = "pluritest_image$num.png";
                my $plot_path = $self->output_file(output_key => 'pluritest_plots', basename => $basename, type => 'bin', metadata => $meta)->path->stringify;
                $self->relate_input_to_output($conv_path, 'pluritest_plot', $plot_path);
                my $plot_node = $vrtrack->get_file($plot_path);
                $plot_node->add_properties({ type => $num_to_type{$num} });
                $donor_node->relate_to($plot_node, 'pluritest_plot');
            }
            
            my $cmd    = $self->r_cmd_prefix . qq[ --slave --args $conv_path $pluritest_data < $pluritest_script];
            my $vr_cmd = "use VRPipe::Steps::pluritest_plot_gene_expression; VRPipe::Steps::pluritest_plot_gene_expression->plot(cmd_line => q[$cmd], csv_out_path => q[$csv_path]);";
            $self->dispatch_vrpipecode($vr_cmd, $req);
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
    
    method plot (ClassName|Object $self: Str :$cmd_line, Str|File :$csv_out_path!) {
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $csv_out_path);
        $output_file->update_stats_from_disc;
        my $vrtrack              = VRPipe::Schema->create('VRTrack');
        my $output_file_in_graph = $vrtrack->get_file($csv_out_path);
        
        # parse to store per-sample results in nodes attached to each sample
        my $md5 = $output_file_in_graph->md5;
        unless ($md5) {
            $md5 = $output_file->file_md5($output_file);
            $output_file_in_graph->md5($md5);
        }
        
        my $graph = $vrtrack->graph;
        my $sample_source;
        
        my $fh     = $output_file->openr;
        my $header = <$fh>;
        chomp($header);
        my @headings;
        foreach my $col (split(/,/, $header)) {
            $col =~ s/^"|"$//g;
            push(@headings, $col);
        }
        
        my $t = time();
        while (<$fh>) {
            chomp;
            my @cols = split(/,/, $_);
            my $sample = $cols[0];
            $sample =~ s/^"|"$//g;
            $sample_source ||= $vrtrack->sample_source($sample);
            my $sample_props = $vrtrack->sample_props_from_string($sample, $sample_source);
            
            my $data = {};
            foreach my $i (1 .. $#cols) {
                $data->{ $headings[$i] } = $cols[$i];
            }
            
            $vrtrack->add(
                'Pluritest',
                {
                    md5_sample => $md5 . '_' . $sample_props->{name},
                    date       => $t,
                    data       => $graph->json_encode($data)
                },
                incoming => [
                    { node_spec => { namespace => 'VRTrack', label => 'Sample', properties => $sample_props }, type => 'pluritest' },
                    { node => $output_file_in_graph, type => 'parsed' }
                ],
                enqueue => 1
            );
        }
        $output_file->close();
        $vrtrack->dispatch_queue;
    }
}

1;
