use VRPipe::Base;

class VRPipe::Steps::fastq_map with VRPipe::StepRole {
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq', max_files => 2, description => '1-2 fastq files') };
    }
    method body_sub {
        return sub { my $self = shift; };
    }
    method outputs_definition {
        return { bam_file => VRPipe::StepIODefinition->get(type => 'fq', max_files => 1, description => 'mapped bam file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Maps the input fastq(s) with a mapper based on the sequencing technology";
    }
}

1;
