use VRPipe::Base;

class VRPipe::Pipelines::linked_pipeline with VRPipe::PipelineRole {
    method name {
        return 'linked_pipeline';
    }
    
    method description {
        return 'a test pipeline for testing vrpipe sources';
    }
    
    method step_names {
        (
            'text_merge',
            'test_step_one'
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'input_text_files' },
            { from_step => 1, to_step => 2, from_key => 'merged_file', to_key => 'one_input' }
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup',       default_regulation => 1 },
            { after_step => 1, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'delete_inputs', default_regulation => 0 }
        );
    }
}

1;
