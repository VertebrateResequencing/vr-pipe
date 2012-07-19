use VRPipe::Base;

class VRPipe::Pipelines::convex_cnv_calling with VRPipe::PipelineRole {
    method name {
        return 'convex_cnv_calling';
    }
    
    method _num_steps {
        return 2;
    }
    
    method description {
        return 'Run CoNVex pipeline to Generate CNV calls from Read Depth and L2R files';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
             VRPipe::Step->get(name => 'convex_gam_correction'), #
             VRPipe::Step->get(name => 'convex_cnv_call'),       #
            ],
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'rd_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'gam_files', to_key => 'gam_files'),],
            [],);
    }
}

1;
