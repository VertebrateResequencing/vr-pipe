
=head1 NAME

VRPipe::Steps::breakdancer_sv_detection - a step

=head1 DESCRIPTION

Runs the Breakdancer SV detection program using a pre-generated config file

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

class VRPipe::Steps::breakdancer_sv_detection with VRPipe::StepRole {
    method options_definition {
        return { breakdancer_max_options => VRPipe::StepOption->create(description => 'breakdancer_max options excluding bam config file name'),
                 breakdancer_max_exe     => VRPipe::StepOption->create(description => 'full path to breakdancer_max executable', optional => 1, default_value => 'breakdancer_max'),
                 whole_genome_mode       => VRPipe::StepOption->create(description => "Indicates process run in one job, set to 0 to split into seperate jobs by chromosome", optional => 1, default_value => 1),
                 chrom_list              => VRPipe::StepOption->create(description => 'Names of chromosomes if running seperate jobs per chromosome', optional => 1, default_value => '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y'), };
    }
    
    method inputs_definition {
        return { bam_cfg => VRPipe::StepIODefinition->create(type        => 'txt',
                                                             description => 'breakdancer bam config files',
                                                             max_files   => -1) };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options                 = $self->options;
            my $breakdancer_max_exe     = $options->{breakdancer_max_exe};
            my $breakdancer_max_options = $options->{'breakdancer_max_options'};
            my $whole_genome_mode       = $options->{whole_genome_mode};
            my $chrom_list              = $options->{chrom_list};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $bam_cfg (@{ $self->inputs->{bam_cfg} }) {
                my $input_path = $bam_cfg->path;
                my $basename   = $bam_cfg->basename;
                $basename =~ s/\.cfg$//;
                
                if ($whole_genome_mode) {
                    my $breakdancer_max = $self->output_file(output_key => 'breakdancer_max', basename => "${basename}.max", type => 'txt');
                    my $output_path = $breakdancer_max->path;
                    
                    my $cmd = "$breakdancer_max_exe $breakdancer_max_options $input_path > $output_path";
                    $self->dispatch_wrapped_cmd('VRPipe::Steps::breakdancer_sv_detection', 'run_breakdancer_max', [$cmd, $req, { output_files => [$breakdancer_max] }]);
                }
                else {
                    for my $chr (split(' ', $chrom_list)) {
                        my $breakdancer_max = $self->output_file(output_key => 'breakdancer_max', basename => "${basename}.${chr}.max", type => 'txt');
                        my $output_path = $breakdancer_max->path;
                        
                        my $cmd = "$breakdancer_max_exe $breakdancer_max_options -o $chr $input_path > $output_path";
                        $self->dispatch_wrapped_cmd('VRPipe::Steps::breakdancer_sv_detection', 'run_breakdancer_max', [$cmd, $req, { output_files => [$breakdancer_max] }]);
                    }
                }
            }
        };
    }
    
    method outputs_definition {
        return { breakdancer_max => VRPipe::StepIODefinition->create(type        => 'txt',
                                                                     description => 'breakdancer max sv detection results',
                                                                     max_files   => -1) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Generates breakdancer sv detection results from input bam config file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method run_breakdancer_max (ClassName|Object $self: Str $cmd_line) {
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my ($output_path) = $cmd_line =~ /.* > (\S+)$/;
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        
        if ($output_file->num_records == 0) {
            $output_file->unlink;
            $self->throw("Output $output_path is empty)");
        }
        else {
            return 1;
        }
    }

}

1;
