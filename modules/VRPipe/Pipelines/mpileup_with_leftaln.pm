use VRPipe::Base;

class VRPipe::Pipelines::mpileup_with_leftaln with VRPipe::PipelineRole {
    method name {
        return 'mpileup_with_leftaln';
    }
    method _num_steps {
        return 2;
    }
    method description {
        return 'Run samtools mpileup to generate vcf files and left_align indels';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return (
            [ 
                VRPipe::Step->get(name => 'mpileup_vcf'), #1
                VRPipe::Step->get(name => 'vcf_index'),   #2
                VRPipe::Step->get(name => 'gatk_vcf_leftalign'), #3
            ],
   	        [
                VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2,  from_key => 'vcf_files', to_key => 'vcf_files'),
                VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3,  from_key => 'vcf_files', to_key => 'vcf_files'),
            ],
            [ 
                VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1, 2], regulated_by => 'cleanup', default_regulation => 1),
            ]
        );
    }
}

1;
