use VRPipe::Base;

class VRPipe::Pipelines::convex_l2r_bp_generation with VRPipe::PipelineRole {
    method name {
        return 'convex_l2r_bp_generation';
    }
    method _num_steps {
        return 2;
    }
    method description {
        return 'Run CoNVex pipeline to Generate L2R files from Read Depth files, and Breakpoint file, for subsequent CNV Calling pipeline';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ 
					VRPipe::Step->get(name => 'convex_breakpoints'),   #
					VRPipe::Step->get(name => 'convex_L2R'),           #
				],
   	            [ 
					VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'rd_files'),
					VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'rd_files'),
   	            ],
   	            [ ],
				);
    }
}

1;
