use VRPipe::Base;

class VRPipe::Pipelines::vcf_vep_annotate with VRPipe::PipelineRole {
    method name {
        return 'vcf_vep_annotate';
    }
    method _num_steps {
        return 2;
    }
    method description {
        return 'Annotate VCF files with consequences using VEP';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'vep_analysis'), 
                  VRPipe::Step->get(name => 'vcf_vep_consequences') ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'vcf_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'vep_txt', to_key => 'vep_txt'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'vcf_files') ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;
