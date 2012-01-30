use VRPipe::Base;

class VRPipe::Pipelines::genotype_checking_exome with VRPipe::PipelineRole {
    method name {
        return 'genotype_checking_exome';
    }
    method _num_steps {
        return 4;
    }
    method description {
        return 'Performs genotype checking for bam files of exome samples.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'bam_snp_sites'),#1
                  VRPipe::Step->get(name => 'mpileup_bcf_snp_sites'),#2
                  VRPipe::Step->get(name => 'glf_check_genotype'),#3
                  VRPipe::Step->get(name => 'gtypex_genotype_analysis'),#4
                  ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'bam_files_with_metadata', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'bcf_files_with_metadata', to_key => 'bcf_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'gtypex_files_with_metadata', to_key => 'gtypex_files'),
                   ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 4, behaviour => 'delete_outputs', act_on_steps => [1,2], regulated_by => 'cleanup', default_regulation => 0),
                 ]);
    }
}

1;