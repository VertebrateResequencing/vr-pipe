use VRPipe::Base;

class VRPipe::Pipelines::convex_read_depth_generation with VRPipe::PipelineRole {
    method name {
        return 'convex_read_depth_generation';
    }
    
    method _num_steps {
        return 2;
    }
    
    method description {
        return 'Run CoNVex pipeline to Generate Read Depth files and add Sample metadata to bam files, for subsequent L2R Generation and CNV Calling pipelines';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
             VRPipe::Step->get(name => 'convex_read_depth'),     #
             VRPipe::Step->get(name => 'bam_metadata_with_sex'), #
            ],
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'),],
            [],);
    }
}

1;
