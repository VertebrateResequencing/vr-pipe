use VRPipe::Base;

class VRPipe::Pipelines::vcf_chunked_vep_annotate with VRPipe::PipelineRole {
	method name {
		return 'vcf_chunked_vep_annotate';
	}
	method _num_steps {
		return 7;
	}
	method description {
		return 'Run vcf consequence annotation using VEP with vcf file chunking';
	}
	method steps {
		$self->throw("steps cannot be called on this non-persistent object");
	}

	method _step_list {
		return ([
				VRPipe::Step->get(name => 'chunk_genomic_region'), #1
				VRPipe::Step->get(name => 'vcf_split'),            #2
				VRPipe::Step->get(name => 'vcf_annotate'),         #3
				VRPipe::Step->get(name => 'vep_analysis'),         #4
				VRPipe::Step->get(name => 'vcf_vep_consequences'), #5
				VRPipe::Step->get(name => 'vcf_concat'),           #6
				VRPipe::Step->get(name => 'vcf_index'),            #7
				],
				[ 
				VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'vcf_files'),
				VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'chunked_regions_file', to_key => 'chunked_regions_file'), 
				VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'vcf_split_files', to_key => 'vcf_files'),
				VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'annotated_vcf', to_key => 'vcf_files'),
				VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 5, from_key => 'annotated_vcf', to_key => 'vcf_files'),
				VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'vep_txt', to_key => 'vep_txt'),
				VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'conseq_vcf', to_key => 'vcf_files'), 
				VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'concat_vcf', to_key => 'vcf_files'),
				],
				[ VRPipe::StepBehaviourDefiner->new(after_step => 6, behaviour => 'delete_outputs', act_on_steps => [1,2,3,4,5,], regulated_by => 'cleanup', default_regulation => 1) ]	
			   );
	}
}

1;
