use VRPipe::Base;

class VRPipe::Pipelines::gatk_variant_calling_and_filter_vcf with VRPipe::PipelineRole {
    method name {
        return 'gatk_variant_calling_and_filter_vcf';
    }
    method _num_steps {
        return 4;
    }
    method description {
        return 'Call variants with the GATK universal genotyper, then hard-filter the results with GATK variant filtration';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'bam_index'),#1
		  VRPipe::Step->get(name => 'gatk_genotype'),#2
		  VRPipe::Step->get(name => 'vcf_index'),#3
		  VRPipe::Step->get(name => 'gatk_variant_filter'),#4
		  VRPipe::Step->get(name => 'vcf_index')#5
		],
		
		[ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
		  VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'),
		  VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'vcf_files', to_key => 'vcf_files'),
		  VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 4, from_key => 'vcf_files', to_key => 'vcf_files'),
		  VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'filtered_vcf_files', to_key => 'vcf_files')
		],
		
		[ VRPipe::StepBehaviourDefiner->new(after_step => 4, behaviour => 'delete_outputs', act_on_steps => [2, 3], regulated_by => 'delete_unfiltered_vcfs', default_regulation => 1) ]);
    }
}

1;
