
=head1 NAME

VRPipe::Steps::kallisto_index - a step

=head1 DESCRIPTION

Builds a kallisto index of a transcriptome.

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2016 Genome Research Limited.

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

class VRPipe::Steps::kallisto_index with VRPipe::StepRole {
    method options_definition {
        return {
            kallisto_exe => VRPipe::StepOption->create(
                description   => 'path to your kallisto executable',
                optional      => 1,
                default_value => 'kallisto'
            ),
            kallisto_index_options => VRPipe::StepOption->create(
                description   => 'options to kallisto index, excluding -i',
                optional      => 1,
                default_value => '--kmer-size=31'
            ),
            transcripts_fasta => VRPipe::StepOption->create(description => 'absolute path to transcriptome file in FASTA format'),
        };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $kallisto_exe  = $options->{kallisto_exe};
            my $kallisto_opts = $options->{kallisto_index_options};
            if ($kallisto_opts =~ /\s-i/) {
                $self->throw("kallisto_index_options should not include --index=");
            }
            my $fasta = file($options->{transcripts_fasta});
            $self->throw("transcripts_fasta must be an absolute path, $fasta") unless $fasta->is_absolute;
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => "$kallisto_exe",
                    version => VRPipe::StepCmdSummary->determine_version($kallisto_exe, '^kallisto (.+)$'),
                    summary => "kallisto index $kallisto_opts -i \$fasta_index \$fasta"
                )
            );
            
            $self->output_file(output_key => 'transcriptome_index', output_dir => $fasta->dir, basename => $fasta->basename . "_kallisto_index", type => 'bin');
            my $req = $self->new_requirements(memory => 5000, time => 1);
            $self->dispatch(["cd " . $fasta->dir . "; $kallisto_exe index $kallisto_opts -i $fasta" . "_kallisto_index $fasta", $req, { block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return {
            transcriptome_index => VRPipe::StepIODefinition->create(type => 'bin', description => 'transcriptome index file'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Builds an index from a FASTA formatted file of target transcriptome sequences.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
