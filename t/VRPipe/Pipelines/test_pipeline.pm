use VRPipe::Base;

class VRPipe::Pipelines::test_pipeline with VRPipe::PipelineRole {
    method name {
        return 'test_pipeline';
    }
    method _num_steps {
        return 4;
    }
    method description {
        return 'a test pipeline for testing';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'test_step_one'),
                  VRPipe::Step->get(name => 'test_step_two'),
                  VRPipe::Step->get(name => 'test_step_three'),
                  VRPipe::Step->get(name => 'test_step_four'),
                 ],
                 [ {before_step_number => 1, adaptor_hash => { one_input => 'data_element' }},
                   {before_step_number => 2, adaptor_hash => { two_input => 'one_output' }},
                   {before_step_number => 3, adaptor_hash => { three_input => 'two_output' }},
                   {before_step_number => 4, adaptor_hash => { four_input => 'three_output' }} ],
                 [ ]);
    }
}

1;