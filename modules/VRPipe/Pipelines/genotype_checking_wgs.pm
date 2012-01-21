use VRPipe::Base;

class VRPipe::Pipelines::genotype_checking_wgs with VRPipe::PipelineRole {
    method name {
        return 'genotype_checking_wgs';
    }
    method _num_steps {
        return 3;
    }
    method description {
        return 'Performs genotype checking for bam files of whole genome study samples.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'genotype_mpileup_wgs'),#1
                  VRPipe::Step->get(name => 'genotype_checking'),#2
                  VRPipe::Step->get(name => 'genotype_analysis'),#3
                  ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'bcf_files_with_metadata', to_key => 'bcf_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'gtypex_files_with_metadata', to_key => 'gtypex_files'),
                   ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 3, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1),]);
    }
}

1;