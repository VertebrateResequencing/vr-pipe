
=head1 NAME

VRPipe::Steps::plot_polysomy - a step

=head1 DESCRIPTION

This step runs plot-polysomy to generate per-sample copy numbers for each 
chromosome.

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Steps::plot_polysomy with VRPipe::StepRole {
    use VRPipe::Schema;
    
    method options_definition {
        return {
            plot_polysomy_exe => VRPipe::StepOption->create(
                description   => 'path to your plot-polysomy script',
                optional      => 1,
                default_value => 'plot-polysomy.py'
            ),
        };
    }
    
    method inputs_definition {
        return {
            dist_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => '1 or more dist.dat files from polysomy',
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self              = shift;
            my $options           = $self->options;
            my $plot_polysomy_exe = $options->{plot_polysomy_exe};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'plot-polysomy',
                    version => 0,
                    summary => "plot-polysomy -o \$png_file [title\@dir ...]"
                )
            );
            
            my $merged_meta = $self->combined_metadata($self->inputs->{dist_files});
            my $png_file_path = $self->output_file(output_key => 'png_file', basename => 'copy_numbers.png', type => 'png', metadata => $merged_meta)->path->stringify;
            
            my (@args, $sample_source);
            my $vrtrack              = VRPipe::Schema->create('VRTrack');
            my $output_file_in_graph = $vrtrack->add_file($png_file_path);
            foreach my $file (@{ $self->inputs->{dist_files} }) {
                my $dist  = $file->metadata;
                my $query = $dist->{sample};
                push(@args, "$query\@" . $file->dir);
                
                $self->relate_input_to_output($file->path->stringify, 'copy_number_plot', $png_file_path);
                
                # also attach it to the sample, deleting any existing
                # copy_number_by_chromosome_plot relationship first
                $sample_source ||= $vrtrack->sample_source($query);
                my $sample_props = $vrtrack->sample_props_from_string($query, $sample_source);
                my $sample_node = $vrtrack->get('Sample', $sample_props);
                unless ($sample_node) {
                    $self->throw("No sample node with name $sample_props->{name} was found in the graph db (for file " . $file->path . ", query $query)");
                }
                my @existing_plots = $sample_node->related(outgoing => { type => 'copy_number_by_chromosome_plot' });
                foreach my $plot (@existing_plots) {
                    $sample_node->divorce_from($plot);
                }
                $sample_node->relate_to($output_file_in_graph, 'copy_number_by_chromosome_plot');
            }
            
            my $req = $self->new_requirements(memory => 1000, time => 1);
            my $cmd_line = "$plot_polysomy_exe -o $png_file_path @args";
            $self->dispatch([$cmd_line, $req]);
        };
    }
    
    method outputs_definition {
        return {
            png_file => VRPipe::StepIODefinition->create(
                type        => 'png',
                min_files   => 1,
                max_files   => 1,
                description => 'output copy number plot from polysomy',
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run plot-polysomy to generate copy number plots by chromosome.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }

}

1;
