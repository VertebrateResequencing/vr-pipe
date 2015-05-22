
=head1 NAME

VRPipe::Steps::fmd_index - a step

=head1 DESCRIPTION

Build an FMD index for a set of reads using ropebwt2

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

class VRPipe::Steps::fmd_index with VRPipe::StepRole {
    method options_definition {
        return {
            ropebwt2_exe     => VRPipe::StepOption->create(description => 'path to your ropebwt2 executable', optional => 1, default_value => 'ropebwt2'),
            ropebwt2_options => VRPipe::StepOption->create(description => 'options to ropebwt2',              optional => 1, default_value => '-dNCr'),
        };
    }
    
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->create(type => 'fq', max_files => -1, description => 'sequence files to be indexed') };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $ropebwt2_exe  = $options->{ropebwt2_exe};
            my $ropebwt2_opts = $options->{ropebwt2_options};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'ropebwt2',
                    version => VRPipe::StepCmdSummary->determine_version($ropebwt2_exe, '^Usage:\s+ropebwt2-(\S+)'),
                    summary => 'ropebwt2 ' . $ropebwt2_opts . ' $reads_file > $fmd_index'
                )
            );
            
            my $req = $self->new_requirements(memory => 35000, time => 1);
            foreach my $fq (@{ $self->inputs->{fastq_files} }) {
                my $prefix = $fq->basename;
                $prefix =~ s/\.(fq|fastq)(\.gz)?//;
                my $fmd_index = $self->output_file(output_key => 'fmd_index_files', basename => "$prefix.fmd", type => 'bin', metadata => $fq->metadata);
                my $this_cmd = "$ropebwt2_exe $ropebwt2_opts " . $fq->path . " > " . $fmd_index->path;
                $self->dispatch([$this_cmd, $req, { output_files => [$fmd_index] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            fmd_index_files => VRPipe::StepIODefinition->create(type => 'bin', description => 'the ropebwt2 FMD index', max_files => -1),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Build an FMD index for a set of reads using ropebwt2";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
