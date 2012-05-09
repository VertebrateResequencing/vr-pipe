use VRPipe::Base;

class VRPipe::Pipelines::vrtrack_auto_qc with VRPipe::PipelineRole {
    method name {
        return 'vrtrack_auto_qc';
    }
    method _num_steps {
        return 1;
    }
    method description {
        return 'Considering the stats in the bamcheck file for a lane, and the metadata stored on the bam file and in the VRTrack database for the corresponding lane, automatically decide if the lane passes the quality check.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'vrtrack_auto_qc') ],
                [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                  VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bamcheck_files') ],
                [ ]);
    }
}

1;