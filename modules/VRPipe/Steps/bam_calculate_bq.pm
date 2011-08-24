use VRPipe::Base;

class VRPipe::Steps::bam_calculate_bq with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->get(description => 'path to your samtools executable',
                                                         optional => 1,
                                                         default_value => 'samtools') };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            $self->throw("bam_calculate_bq not yet implemented");
        };
    }
    method outputs_definition {
        return { bq_bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => 'a bam file with BQ tag and good NM & MD tags') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Corrects NM & MD tags and calculates BQ (BAQ scores) to aid downstream variant calling";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
