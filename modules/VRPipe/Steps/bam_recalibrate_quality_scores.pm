use VRPipe::Base;

class VRPipe::Steps::bam_recalibrate_quality_scores with VRPipe::StepRole {
    method options_definition {
        return { };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more coordinate-sorted bam files') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            $self->throw("bam_recalibrate_quality_scores not yet implemented");
        };
    }
    method outputs_definition {
        return { recalibrated_bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => 'a bam file with recalibrated quality scores; OQ tag holds the original quality scores') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Recalibrate quality scores of each mapped base using GATK";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
