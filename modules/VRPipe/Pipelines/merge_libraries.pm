use VRPipe::Base;

class VRPipe::Pipelines::merge_libraries with VRPipe::PipelineRole {
    method name {
        return 'merge_libraries';
    }
    method _num_steps {
        return 3;
    }
    method description {
        return 'Merge and reheader';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'bam_merge'),#1
                  VRPipe::Step->get(name => 'sequence_dictionary'),#2
                  VRPipe::Step->get(name => 'bam_reheader'),#3
                  ],

                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3, from_key => 'merged_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'reference_dict', to_key => 'dict_file'),
                  ],

                 [ 
                   VRPipe::StepBehaviourDefiner->new(after_step => 3, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 3, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'cleanup', default_regulation => 0), # should this be regulated by a different option?
                 ]);
    }
}

1;