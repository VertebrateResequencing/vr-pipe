use VRPipe::Base;

class VRPipe::Pipelines::snp_calling_mpileup_bcf with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_mpileup_bcf';
    }
    method _num_steps {
        return 2;
    }
    method description {
        return 'Run samtools mpileup generating both bcf and vcf files';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
		return ([
				VRPipe::Step->get(name => 'mpileup_bcf'),
				VRPipe::Step->get(name => 'bcf_to_vcf'),
				],
				[
				VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
				VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'bcf_files', to_key => 'bcf_files'),
				],
				[
				VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1),
				]);
    }
}

1;
