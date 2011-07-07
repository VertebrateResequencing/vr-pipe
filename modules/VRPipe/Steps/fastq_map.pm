use VRPipe::Base;

class VRPipe::Steps::fastq_map with VRPipe::StepRole {
    method options_definition {
        return { };
    }
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq', max_files => 3, description => '1-3 fastq files') };
    }
    method body_sub {
        return sub { my $self = shift; }; #*** wrappers return command lanes that you then dispatch as normal
    }
    method outputs_definition {
        return { mapped_bam_files => VRPipe::StepIODefinition->get(type => 'fq', max_files => 2, description => 'mapped bam file(s)') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Maps the input fastq(s) with a mapper based on the sequencing technology";
    }
}

1;
