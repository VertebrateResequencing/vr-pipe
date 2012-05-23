use VRPipe::Base;

class VRPipe::Pipelines::bam_index with VRPipe::PipelineRole {
    method name {
        return 'bam_index';
    }
    method _num_steps {
        return 1;
    }
    method description {
        return 'Index bam files';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'bam_index') ],
                [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files') ],
                [ ]);
    }
}

1;