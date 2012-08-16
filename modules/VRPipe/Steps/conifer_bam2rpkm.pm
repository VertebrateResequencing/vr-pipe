
=head1 NAME

VRPipe::Steps::conifer_bam2rpkm - a step

=head1 DESCRIPTION

Generates conifer RPKM (reads per thousand bases per million reads) text files
from a set of input bams, to be analysed by conifer as a group

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

class VRPipe::Steps::conifer_bam2rpkm with VRPipe::StepRole {
    method options_definition {
        return {
            python_exe  => VRPipe::StepOption->create(description => 'full path to python executable', optional => 1, default_value => 'python'),
            conifer_py  => VRPipe::StepOption->create(description => 'full path to conifer.py',        optional => 1, default_value => 'conifer.py'),
            probes_file => VRPipe::StepOption->create(description => 'probes / target definition file'),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                description => 'set of bam files to be analysed as a group',
                max_files   => -1
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options     = $self->options;
            my $python_exe  = $options->{python_exe};
            my $conifer_py  = $options->{conifer_py};
            my $probes_file = $options->{probes_file};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam_file (@{ $self->inputs->{bam_files} }) {
                my $basename = $bam_file->basename;
                $basename =~ s/\.bam$/.rpkm.txt/;
                my $rpkm_out = $self->output_file(output_key => 'rpkm_out', basename => $basename, type => 'txt');
                
                my $input_path  = $bam_file->path;
                my $output_path = $rpkm_out->path;
                
                my $cmd = "$python_exe $conifer_py rpkm --probes $probes_file --input $input_path --output $output_path";
                $self->warn($cmd);
                $self->dispatch_wrapped_cmd('VRPipe::Steps::conifer_bam2rpkm', 'gen_rpkm', [$cmd, $req, { output_files => [$rpkm_out] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            rpkm_out => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'one conifer rpkm output file per bam',
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Generates conifer RPMK (reads per thousand bases per million reads) text files from a set of input bams to be analysed by conifer as a group";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method gen_rpkm (ClassName|Object $self: Str $cmd_line) {
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my ($output_path) = $cmd_line =~ /--output (\S+)$/;
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        
        if ($output_file->lines == 0) {
            $output_file->unlink;
        }
        else {
            return 1;
        }
    }
}

1;
