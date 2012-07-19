use VRPipe::Base;

class VRPipe::Pipelines::linked_pipeline with VRPipe::PipelineRole {
    method name {
        return 'linked_pipeline';
    }
    
    method _num_steps {
        return 2;
    }
    
    method description {
        return 'a test pipeline for testing vrpipe sources';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([VRPipe::Step->get(name => 'text_merge'), VRPipe::Step->get(name => 'test_step_one')], [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'input_text_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'merged_file', to_key => 'one_input')], [VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1), VRPipe::StepBehaviourDefiner->new(after_step => 1, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_inputs', default_regulation => 0)]);
    }
}

1;
