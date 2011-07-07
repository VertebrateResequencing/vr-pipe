use VRPipe::Base;

class VRPipe::Steps::bam_stats with VRPipe::StepRole {
    method options_definition {
        return { };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files to calculate stats for') };
    }
    method body_sub {
        return sub { my $self = shift; };
    }
    method outputs_definition {
        return { bas_files => VRPipe::StepIODefinition->get(type => 'fq', max_files => -1, description => 'a .bas file for each input bam file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Calculates various statistics about bam files, producing .bas files";
    }
}

1;
