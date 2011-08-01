use VRPipe::Base;

class VRPipe::Steps::bam_realignment_around_known_indels with VRPipe::StepRole {
    method options_definition {
        return { };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            $self->throw("bam_realignment_around_known_indels not yet implemented");
        };
    }
    method outputs_definition {
        return { realigned_bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => 'a name-sorted uncompressed bam file with improved alignements near indels') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Realigns reads around known indels to improve subsequent variant calling, producing a name-sorted uncompressed bam";
    }
}

1;
