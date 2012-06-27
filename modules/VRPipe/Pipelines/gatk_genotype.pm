use VRPipe::Base;

class VRPipe::Pipelines::gatk_genotype with VRPipe::PipelineRole {
    method name {
        return 'gatk_genotype';
    }
    method _num_steps {
        return 2;
    }
    method description {
        return 'Run gatk universal genotyper';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ 
					VRPipe::Step->get(name => 'gatk_genotype'),                  #1
					VRPipe::Step->get(name => 'vcf_index'),                      #2
				],
   	            [ 
					VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
					VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'vcf_files', to_key => 'vcf_files'),
				],
   	            [ ],
				);
    }
}

1;
