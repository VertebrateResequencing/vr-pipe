use VRPipe::Base;

class VRPipe::Pipelines::archive_files with VRPipe::PipelineRole {
    method name {
        return 'archive_files';
    }
    method _num_steps {
        return 1;
    }
    method description {
        return 'Safely move files from one disc to a pool of one or more other discs, eg. for archival purposes. (the output_root option for this pipeline is meaningless)';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'archive_files') ],
                [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'file') ],
                [ ]);
    }
}

1;