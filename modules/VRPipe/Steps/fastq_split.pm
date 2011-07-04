use VRPipe::Base;

class VRPipe::Steps::fastq_split with VRPipe::StepRole {
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq', max_files => 3, description => '1-3 fastq files') };
    }
    method body_sub {
        return sub { my $self = shift; };
    }
    method outputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq', max_files => -1, description => 'split fastq files') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Takes a single-ended fastq file and/or 2 fastq files of a paired-end read and splits them into multiple smaller fastq files";
    }
}

1;
