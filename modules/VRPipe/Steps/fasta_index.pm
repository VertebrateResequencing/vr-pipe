use VRPipe::Base;

class VRPipe::Steps::fasta_index with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->get(description => 'path to your samtools executable', optional => 1, default_value => 'samtools'),
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file') };
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $samtools = $options->{samtools_exe};
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path, $ref") unless $ref->is_absolute;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'samtools', version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'), summary => 'samtools faidx $reference_fasta'));
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            $self->output_file(output_key => 'fai_file',
                               output_dir => $ref->dir,
                               basename => $ref->basename.'.fai',
                               type => 'txt');
            $self->dispatch([qq[$samtools faidx $ref], $req, {block_and_skip_if_ok => 1}]);
        };
    }
    method outputs_definition {
        return { fai_file => VRPipe::StepIODefinition->get(type => 'txt', description => 'a .fai index file for the reference fasta') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Indexes fasta files using samtools";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
