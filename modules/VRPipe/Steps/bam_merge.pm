use VRPipe::Base;

class VRPipe::Steps::bam_merge with VRPipe::StepRole {
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files to merge') };
    }
    method body_sub {
        return sub { my $self = shift; };
    }
    method outputs_definition {
        return { bam_file => VRPipe::StepIODefinition->get(type => 'fq', max_files => 1, description => 'merged bam file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Merges multiple bam files assumed to be for the same lane into a single bam file";
    }
}

1;
