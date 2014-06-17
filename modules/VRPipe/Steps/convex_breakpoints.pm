
=head1 NAME

VRPipe::Steps::convex_breakpoints - a step

=head1 DESCRIPTION

Runs the BreakpointsCall R script from the CoNVex packages. Runs once per
pipeline, generating a Breakpoints file from a sample read depth file.

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

class VRPipe::Steps::convex_breakpoints extends VRPipe::Steps::r_script {
    around options_definition {
        return {
            %{ $self->$orig },
            'rd_sample_file_name' => VRPipe::StepOption->create(description => 'Full path to a specific sample Read Depth file from which to generate breakpoints', optional => 1),
            'convex_rscript_path' => VRPipe::StepOption->create(description => 'Full path to CoNVex R scripts'),
            'max_bin_size' => VRPipe::StepOption->create(description => 'Maximum bin size', optional => 1, default_value => 1000),
        };
    }
    
    method inputs_definition {
        return { rd_files => VRPipe::StepIODefinition->create(type => 'rd', max_files => -1, description => 'Set of convex Read Depths files (eg grouped by metadata)', metadata => { sample => 'sample name' }), };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $rd_sample_file_name = $options->{'rd_sample_file_name'};
            my $convex_rscript_path = $options->{'convex_rscript_path'};
            my $max_bin_size        = $options->{'max_bin_size'};
            
            my $req = $self->new_requirements(memory => 1000, time => 1);
            
            # Select first Read Depth file in pipeline if none provided as a param
            my $sample_rd_file;
            if ($rd_sample_file_name) {
                my $rd_file = file($rd_sample_file_name);
                $self->throw("rd_sample_file_name must be absolute path") unless $rd_file->is_absolute;
                $sample_rd_file = VRPipe::File->create(path => $rd_sample_file_name);
            }
            else {
                $sample_rd_file      = $self->inputs->{rd_files}[0];
                $rd_sample_file_name = $sample_rd_file->path;
            }
            $self->throw(" not read depth sample file chosen") unless $sample_rd_file;
            
            my $breakpoints_file = $self->output_file(output_key => 'breakpoints_file', output_dir => $sample_rd_file->dir, basename => "breakpoints.bp", type => 'bp', metadata => { source_rd_file => $sample_rd_file->path });
            my $breakpoints_path = $breakpoints_file->path;
            
            my $cmd = $self->rscript_cmd_prefix . " $convex_rscript_path/BreakpointsCall.R $rd_sample_file_name,$max_bin_size,$breakpoints_path";
            
            $self->dispatch([$cmd, $req]);
        };
    }
    
    method outputs_definition {
        return {
            breakpoints_file => VRPipe::StepIODefinition->create(type => 'bp', max_files => 1, description => 'breakpoints file', metadata => { source_rd_file => 'rd file used to create breakpoints file' }),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs the BreakpointsCall R script from the CoNVex packages, generating a Breakpoints file from a sample read depth file";
    }
    
    method max_simultaneous {
        return 1;            # should only run once
    }
}

1;
