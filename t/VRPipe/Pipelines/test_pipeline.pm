use VRPipe::Base;

class VRPipe::Pipelines::test_pipeline with VRPipe::PipelineRole {
    method name {
        return 'test_pipeline';
    }
    
    method description {
        return 'a test pipeline for testing';
    }
    
    method step_names {
        (
            'test_step_one',
            'test_step_two',
            'test_step_three',
            'test_step_four'
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'one_input' },
            { from_step => 1, to_step => 2, from_key => 'one_output', to_key => 'two_input' },
            { from_step => 2, to_step => 3, from_key => 'two_output', to_key => 'three_input' },
            { from_step => 3, to_step => 4, from_key => 'three_output', to_key => 'four_input' }
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 3, behaviour => 'delete_outputs', act_on_steps => [2], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 4, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'cleanup', default_regulation => 1 }
        );
    }
}

1;
