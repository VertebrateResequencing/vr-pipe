use VRPipe::Base;

class VRPipe::Pipelines::conifer_analysis with VRPipe::PipelineRole {
    method name {
        return 'conifer_analysis';
    }
    method _num_steps {
        return 2;
    }
    method description {
        return 'Runs conifer analysis for grouped sets of exome bams, generating a SVD-ZRPKM values, cnv calls and associated plots for each group',
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ 
                  VRPipe::Step->get(name => 'conifer_bam2rpkm'), 
                  VRPipe::Step->get(name => 'conifer_analyze'),
                ],
                 
                [ 
                  VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                  VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'rpkm_out', to_key => 'rpkm_in'),
                ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;
