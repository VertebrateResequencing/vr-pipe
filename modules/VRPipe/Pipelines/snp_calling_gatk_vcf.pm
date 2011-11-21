use VRPipe::Base;

class VRPipe::Pipelines::snp_calling_gatk_vcf with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_gatk_vcf';
    }
    method _num_steps {
        return 1;
    }
    method description {
        return 'Run gatk universal genotyper to generate vcf files';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'gatk_genotype') ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files') ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 1, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1) ]	
				);
    }
}

1;
