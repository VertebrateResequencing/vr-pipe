use VRPipe::Base;

class VRPipe::Pipelines::import_bams with VRPipe::PipelineRole {
    method name {
        return 'import_bams';
    }
    
    method _num_steps {
        return 2;
    }
    
    method description {
        return 'a pipeline for importing bams, getting their metadata';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
                VRPipe::Step->get(name => 'test_import_bams'), #1
                VRPipe::Step->get(name => 'bam_metadata'),     #2
            ],
            
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'local_bam_files', to_key => 'bam_files')],
            
            []
        );
    }
}

1;
