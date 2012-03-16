use VRPipe::Base;

class VRPipe::Pipelines::upload_ega with VRPipe::PipelineRole {
    method name {
        return 'upload_ega';
    }
    method _num_steps {
        return 1;
    }
    method description {
        return 'Uploads files into the EGA dropbox using the EGA upload client application';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'ega_upload') ],
   	        [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'upload_files') ],
                [ ]);
    }
}

1;
