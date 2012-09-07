
=head1 NAME

VRPipe::Steps::fasta_index - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::fasta_index with VRPipe::StepRole {
    method options_definition {
        return {
            samtools_exe    => VRPipe::StepOption->create(description => 'path to your samtools executable', optional => 1, default_value => 'samtools'),
            reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file')
        };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options  = $self->options;
            my $samtools = $options->{samtools_exe};
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path, $ref") unless $ref->is_absolute;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'samtools', version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'), summary => 'samtools faidx $reference_fasta'));
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            $self->output_file(
                output_key => 'fai_file',
                output_dir => $ref->dir,
                basename   => $ref->basename . '.fai',
                type       => 'txt'
            );
            $self->dispatch([qq[$samtools faidx $ref], $req, { block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return { fai_file => VRPipe::StepIODefinition->create(type => 'txt', description => 'a .fai index file for the reference fasta') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Indexes fasta files using samtools";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
