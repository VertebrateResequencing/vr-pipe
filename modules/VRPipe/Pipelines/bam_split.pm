use VRPipe::Base;

class VRPipe::Pipelines::bam_split with VRPipe::PipelineRole {
    method name {
        return 'bam_split';
    }
    method _num_steps {
        return 1;
    }
    method description {
        return 'Split bams into chromosomal bams';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'bam_split_by_sequence')#1
                  ],

                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                  ],

                 [ VRPipe::StepBehaviourDefiner->new(after_step => 1, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0),
                  ]);
    }
}

1;