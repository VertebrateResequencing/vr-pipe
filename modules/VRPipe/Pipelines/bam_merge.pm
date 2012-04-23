use VRPipe::Base;

class VRPipe::Pipelines::bam_merge with VRPipe::PipelineRole {
    method name {
        return 'bam_merge';
    }
    method _num_steps {
        return 2;
    }
    method description {
        return 'Merge bam files together and index the results (suitable for merging library bams to the sample level, and "merge-across")';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'bam_merge'),#1
                  VRPipe::Step->get(name => 'bam_index'),#2
                  ],

                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'merged_bam_files', to_key => 'bam_files'),
                  ],

                 [ VRPipe::StepBehaviourDefiner->new(after_step => 1, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0) ]);
    }
}

1;