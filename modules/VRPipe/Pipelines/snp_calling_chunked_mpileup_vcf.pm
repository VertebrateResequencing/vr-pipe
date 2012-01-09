use VRPipe::Base;

class VRPipe::Pipelines::snp_calling_chunked_mpileup_vcf with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_chunked_mpileup_vcf';
    }
    method _num_steps {
        return 3;
    }
    method description {
        return 'Run samtools mpileup with bam chunking to generate vcf files';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
					VRPipe::Step->get(name => 'chunk_genomic_region'),
					VRPipe::Step->get(name => 'mpileup_vcf'),
					VRPipe::Step->get(name => 'vcf_concat'),
				],
   	            [ 
					VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'),
					VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'chunked_regions_file', to_key => 'chunked_regions_file'), 
					VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'vcf_files', to_key => 'vcf_files'), 
				],
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 3, behaviour => 'delete_outputs', act_on_steps => [2], regulated_by => 'cleanup', default_regulation => 1) ]	
				);
    }
}

1;
