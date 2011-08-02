use VRPipe::Base;

class VRPipe::Steps::bam_fix_mates with VRPipe::StepRole {
    method options_definition {
        return { };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more name-sorted bam files') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            $self->throw("bam_fix_mates not yet implemented");
        };
    }
    method outputs_definition {
        return { fixmate_bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => 'a coordinate-sorted uncompressed bam file with fixed mates') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Fixes mate information and coordinate sorts a name-sorted bam file";
    }
}

1;
