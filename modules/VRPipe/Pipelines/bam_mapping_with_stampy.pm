use VRPipe::Base;

class VRPipe::Pipelines::bam_mapping_with_stampy with VRPipe::PipelineRole {
    method name {
        return 'bam_mapping_with_stampy';
    }
    method _num_steps {
        return 11;
    }
    method description {
        return 'Map reads in bam files to a reference genome with stampy (and bwa)';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'sequence_dictionary'), #1
                  VRPipe::Step->get(name => 'bwa_index'), #2
                  VRPipe::Step->get(name => 'stampy_buildgenome'), #3
                  VRPipe::Step->get(name => 'stampy_buildhash'), #4
                  VRPipe::Step->get(name => 'bam_metadata'), #5
                  VRPipe::Step->get(name => 'bam_name_sort'), #6
                  VRPipe::Step->get(name => 'bam_to_fastq'), #7
                  VRPipe::Step->get(name => 'fastq_split'), #8
                  VRPipe::Step->get(name => 'stampy_map_fastq'), #9
                  VRPipe::Step->get(name => 'sam_to_fixed_bam'), #10
                  VRPipe::Step->get(name => 'bam_merge_lane_splits'), #11
                 ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 5, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 6, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'name_sorted_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 8, from_key => 'fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 8, to_step => 9, from_key => 'split_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 9, to_step => 10, from_key => 'stampy_sam_files', to_key => 'sam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 10, to_step => 11, from_key => 'fixed_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 11, from_key => 'reference_dict', to_key => 'dict_file') ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 7, behaviour => 'delete_outputs', act_on_steps => [6], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 9, behaviour => 'delete_outputs', act_on_steps => [7, 8], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 10, behaviour => 'delete_outputs', act_on_steps => [9], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 11, behaviour => 'delete_outputs', act_on_steps => [10], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;
