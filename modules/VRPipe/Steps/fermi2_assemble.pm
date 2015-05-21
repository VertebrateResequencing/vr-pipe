
=head1 NAME

VRPipe::Steps::fermi2_assemble - a step

=head1 DESCRIPTION

Assemble unitigs from an FMD index of reads using fermi2

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::fermi2_assemble with VRPipe::StepRole {
    method options_definition {
        return {
            fermi2_exe              => VRPipe::StepOption->create(description => 'path to your fermi2 executable', optional => 1, default_value => 'fermi2'),
            fermi2_assemble_options => VRPipe::StepOption->create(description => 'options to fermi2 assemble',     optional => 1, default_value => '-l 71 -m 113 -t 16'),
        };
    }
    
    method inputs_definition {
        return { fmd_index_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'sequence files to be indexed') };
    }
    
    method body_sub {
        return sub {
            my $self        = shift;
            my $options     = $self->options;
            my $fermi2_exe  = $options->{fermi2_exe};
            my $fermi2_opts = $options->{fermi2_assemble_options};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'fermi2',
                    version => VRPipe::StepCmdSummary->determine_version($fermi2_exe, '^Version: (.+)$'),
                    summary => "fermi2 assemble $fermi2_opts \$input.fmd | gzip -1 > \$unitigs.fq.gz"
                )
            );
            
            my ($cpus) = $fermi2_opts =~ m/-t\s*(\d+)/;
            my $req = $self->new_requirements(memory => 30000, time => 1, $cpus ? (cpus => $cpus) : ());
            foreach my $fmd (@{ $self->inputs->{fmd_index_files} }) {
                my $prefix = $fmd->basename;
                $prefix =~ s/\.fmd//;
                my $unitigs = $self->output_file(output_key => 'assembled_unitigs', basename => "$prefix.fq.gz", type => 'fq', metadata => $fmd->metadata);
                my $this_cmd = "$fermi2_exe assemble $fermi2_opts " . $fmd->path . " | gzip -1 > " . $unitigs->path;
                $self->dispatch([$this_cmd, $req, { output_files => [$unitigs] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            assembled_unitigs => VRPipe::StepIODefinition->create(type => 'fq', description => 'fastq file containing the fermi2 assembled unitigs', max_files => -1),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Assemble unitigs from an FMD index of reads using fermi2";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
