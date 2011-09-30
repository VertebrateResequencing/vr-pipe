use VRPipe::Base;

class VRPipe::Pipelines::bam_mapping_with_bwa with VRPipe::PipelineRole {
    method name {
        return 'bam_mapping_with_bwa';
    }
    method _num_steps {
        return 7;
    }
    method description {
        return 'Map reads in bam files to a reference genome with bwa';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'sequence_dictionary'), #1
                  VRPipe::Step->get(name => 'bwa_index'), #2
                  VRPipe::Step->get(name => 'bam_metadata'), #3
                  VRPipe::Step->get(name => 'bwa_aln_bam'), #4
                  VRPipe::Step->get(name => 'bwa_sam_using_bam'), #5
                  VRPipe::Step->get(name => 'sam_to_fixed_bam'), #6
                  VRPipe::Step->get(name => 'bam_reheader'), #7
                  #VRPipe::Step->get(name => 'bam_stats') #8
                 ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'bam_files_with_metadata', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'bwa_sai_files', to_key => 'sai_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 5, from_key => 'bam_files_with_metadata', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'bwa_sam_files', to_key => 'sam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'fixed_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 7, from_key => 'reference_dict', to_key => 'dict_file'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 8, from_key => 'merged_lane_bams', to_key => 'bam_files') ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 6, behaviour => 'delete_outputs', act_on_steps => [4, 5], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 7, behaviour => 'delete_outputs', act_on_steps => [6], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;