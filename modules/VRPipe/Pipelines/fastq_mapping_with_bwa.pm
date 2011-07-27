use VRPipe::Base;

class VRPipe::Pipelines::fastq_mapping_with_bwa with VRPipe::PipelineRole {
    method name {
        return 'fastq_mapping_with_bwa';
    }
    method _num_steps {
        return 9;
    }
    method description {
        return 'Map reads in fastq files to a reference genome with bwa';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'fastq_metadata'),
                  VRPipe::Step->get(name => 'fastq_split'),
                  VRPipe::Step->get(name => 'sequence_dictionary'),
                  VRPipe::Step->get(name => 'bwa_index'),
                  VRPipe::Step->get(name => 'bwa_aln_fastq'),
                  VRPipe::Step->get(name => 'bwa_sam'),
                  VRPipe::Step->get(name => 'sam_to_fixed_bam'),
                  VRPipe::Step->get(name => 'bam_merge_lane_splits'),
                 #VRPipe::Step->get(name => 'bam_stats')
                 ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'fastq_files'), 
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'fastq_files_with_metadata', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 5, from_key => 'split_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 6, from_key => 'split_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'bwa_sai_files', to_key => 'sai_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'bwa_sam_files', to_key => 'sam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 8, from_key => 'fixed_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 8, from_key => 'reference_dict', to_key => 'dict_file'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 8, to_step => 9, from_key => 'merged_lane_bams', to_key => 'bam_files') ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 6, behaviour => 'delete_outputs', act_on_steps => [2, 5], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 7, behaviour => 'delete_outputs', act_on_steps => [6], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 8, behaviour => 'delete_outputs', act_on_steps => [7], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;