
=head1 NAME

VRPipe::Steps::fermi2_simplify - a step

=head1 DESCRIPTION

Simplify a unitig assembly using fermi2

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

class VRPipe::Steps::fermi2_simplify with VRPipe::StepRole {
    method options_definition {
        return {
            fermi2_exe              => VRPipe::StepOption->create(description => 'path to your fermi2 executable', optional => 1, default_value => 'fermi2'),
            fermi2_simplify_options => VRPipe::StepOption->create(description => 'options to fermi2 simplify',     optional => 1, default_value => '-CSo 76 -m 113 -T 71'),
        };
    }
    
    method inputs_definition {
        return { assembled_unitigs => VRPipe::StepIODefinition->create(type => 'fq', max_files => -1, description => 'sequence files to be indexed') };
    }
    
    method body_sub {
        return sub {
            my $self        = shift;
            my $options     = $self->options;
            my $fermi2_exe  = $options->{fermi2_exe};
            my $fermi2_opts = $options->{fermi2_simplify_options};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'fermi2',
                    version => VRPipe::StepCmdSummary->determine_version($fermi2_exe, '^Version: (.+)$'),
                    summary => "fermi2 simplify $fermi2_opts \$unitigs.fq.gz | gzip -1 > \$simple_unitigs.fq.gz"
                )
            );
            
            my $req = $self->new_requirements(memory => 40000, time => 1);
            foreach my $utigs (@{ $self->inputs->{assembled_unitigs} }) {
                my $prefix = $utigs->basename;
                $prefix =~ s/\.(fq|fastq)(\.gz)?//;
                my $simple_unitigs = $self->output_file(output_key => 'simplified_unitigs', basename => "$prefix.fq.gz", type => 'fq', metadata => $utigs->metadata);
                my $this_cmd = "$fermi2_exe simplify $fermi2_opts " . $utigs->path . " | gzip -1 > " . $simple_unitigs->path;
                $self->dispatch([$this_cmd, $req, { output_files => [$simple_unitigs] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            simplified_unitigs => VRPipe::StepIODefinition->create(type => 'fq', description => 'fastq file containing the fermi2 simplified unitigs', max_files => -1),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Simplify a unitig assembly using fermi2";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
