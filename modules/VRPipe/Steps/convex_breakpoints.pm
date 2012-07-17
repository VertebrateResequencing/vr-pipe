=head1 NAME

VRPipe::Steps::convex_breakpoints - a step

=head1 DESCRIPTION

Runs the BreakpointsCall R script from the CoNVex packages. Runs once per pipeline, generating a Breakpoints file from a sample read depth file.

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

class VRPipe::Steps::convex_breakpoints with VRPipe::StepRole {
    method options_definition {
        return { 
                 'rd_sample_file_name' => VRPipe::StepOption->create(description => 'Full path to a specific sample Read Depth file from which to generate breakpoints', optional => 1),
                 'rscript_path' => VRPipe::StepOption->create(description => 'Full path to CoNVex R scripts'),
                 'max_bin_size' => VRPipe::StepOption->create(description => 'Maximum bin size', optional => 1, default_value => 1000),
                 'bp_file_name' => VRPipe::StepOption->create(description => 'Full path to the output Breakpoints file'),
        };
    }
    method inputs_definition {
        return { 
            rd_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'Set of convex Read Depths files (eg grouped by metadata)'),
        };
    }
    method body_sub {
        return sub {
            my $self = shift;

            my $options = $self->options;
            my $rd_sample_file_name = $options->{'rd_sample_file_name'};
            my $rscript_path = $options->{'rscript_path'};
            my $max_bin_size = $options->{'max_bin_size'};
            my $bp_file_name = $options->{'bp_file_name'};

            my $req = $self->new_requirements(memory => 1000, time => 1);

            # Select first Read Depth file in pipeline if none provided as a param
            if ($rd_sample_file_name) {
                my $rd_file = Path::Class::File->new($rd_sample_file_name);
                $self->throw("rd_sample_file_name must be absolute path") unless $rd_file->is_absolute;
            }
            else {
                my $rd_file = $self->inputs->{rd_files}[0];
				$rd_sample_file_name = $rd_file->path;
            }

            my $bp_file = Path::Class::File->new($bp_file_name);
            $self->throw("bp_file_name must be absolute path") unless $bp_file->is_absolute;

            my $cmd = "R --vanilla --slave --args '$rd_sample_file_name,$max_bin_size,$bp_file_name' < $rscript_path/BreakpointsCall.R";
            $self->dispatch([$cmd, $req]);
        };
    }
    method outputs_definition {
        return { };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Runs the BreakpointsCall R script from the CoNVex packages, generating a Breakpoints file from a sample read depth file";
    }
    method max_simultaneous {
        return 1; # should only run once
    }
}

1;
