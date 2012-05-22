use VRPipe::Base;

class VRPipe::Pipelines::breakdancer_analysis with VRPipe::PipelineRole {
    method name {
        return 'breakdancer_analysis';
    }
    method _num_steps {
        return 2;
    }
    method description {
        return 'Run breakdancer structural variant detection';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ 
                  VRPipe::Step->get(name => 'breakdancer_bam2cfg'), 
                  VRPipe::Step->get(name => 'breakdancer_sv_detection'),
                ],
                 
                [ 
                  VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                  VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'bam_cfg', to_key => 'bam_cfg'),
                ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;
