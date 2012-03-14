use VRPipe::Base;

class VRPipe::Pipelines::bam_mapping_with_bwa_via_fastq with VRPipe::PipelineRole {
    method name {
        return 'bam_mapping_with_bwa_via_fastq';
    }
    method _num_steps {
        return 12;
    }
    method description {
        return 'Map reads in bam files to a reference genome with bwa fastq alignment';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'sequence_dictionary'), #1
                  VRPipe::Step->get(name => 'bwa_index'), #2
                  VRPipe::Step->get(name => 'bam_metadata'), #3
                  VRPipe::Step->get(name => 'bam_name_sort'), #4
                  VRPipe::Step->get(name => 'bam_to_fastq'),#5
                  VRPipe::Step->get(name => 'fastq_split'),#6
                  VRPipe::Step->get(name => 'bwa_aln_fastq'),#7
                  VRPipe::Step->get(name => 'bwa_sam'),#8
                  VRPipe::Step->get(name => 'sam_to_fixed_bam'),#9
                  VRPipe::Step->get(name => 'bam_merge_lane_splits'),#10
                  VRPipe::Step->get(name => 'bam_index'),#11
                  VRPipe::Step->get(name => 'bam_stats'),#12
                 ],
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'bam_files_with_metadata', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'name_sorted_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'split_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 8, from_key => 'bwa_sai_files', to_key => 'sai_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 8, from_key => 'split_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 8, to_step => 9, from_key => 'bwa_sam_files', to_key => 'sam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 9, to_step => 10, from_key => 'fixed_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 10, from_key => 'reference_dict', to_key => 'dict_file'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 10, to_step => 11, from_key => 'merged_lane_bams', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 10, to_step => 12, from_key => 'merged_lane_bams', to_key => 'bam_files'),
                 ],
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 8, behaviour => 'delete_outputs', act_on_steps => [4, 5, 6, 7], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 11, behaviour => 'delete_outputs', act_on_steps => [8, 9], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;
