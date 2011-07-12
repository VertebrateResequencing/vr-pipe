use VRPipe::Base;

class VRPipe::Steps::bwa_index with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'genome reference file to map against'),
                 bwa_index_cmd => VRPipe::StepOption->get(description => 'the near-complete bwa index command line, including desired options, but excluding the reference fasta file',
                                                          optional => 1,
                                                          default_value => 'bwa index -a bwtsw')};
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub {
            my $self = shift;
            $self->throw("foo");
        };
    }
    method outputs_definition {
        return { bwa_aln_foo => VRPipe::StepIODefinition->get(type => 'dat', description => 'one of the index files produced by bwa index') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "";
    }
}

1;
