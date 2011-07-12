use VRPipe::Base;

class VRPipe::Steps::bam_merge_lane_splits with VRPipe::StepRole {
    method options_definition {
        return { bam_merge_keep_single_paired_separate => VRPipe::StepOption->get(description => 'when merging bam files, separately merges single ended bam files and paired-end bam files, resulting in 2 merged bam files',
                                                                                  optional => 1,
                                                                                  default_value => 1) };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files to merge') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            $self->throw("foo");
        };
    }
    method outputs_definition {
        return { merged_lane_bams => VRPipe::StepIODefinition->get(type => 'fq', max_files => 2, description => 'merged bam file(s)') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Merges multiple bam files for the same lane into a single bam file (per library layout)";
    }
}

1;
