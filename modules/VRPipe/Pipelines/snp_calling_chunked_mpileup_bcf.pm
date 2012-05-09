use VRPipe::Base;

class VRPipe::Pipelines::snp_calling_chunked_mpileup_bcf with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_chunked_mpileup_bcf';
    }
    method _num_steps {
        return 8;
    }
    method description {
        return 'Create all-sites bcf using samtools mpileup and call variants using bcftools. Accounts for ploidy.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
            return ([ VRPipe::Step->get(name => 'chunk_genomic_region'), #1
                      VRPipe::Step->get(name => 'bam_index'),            #2
                      VRPipe::Step->get(name => 'mpileup_chunked_bcf'),  #3
                      VRPipe::Step->get(name => 'sex_to_ploidy'),        #4
                      VRPipe::Step->get(name => 'bcf_to_vcf'),           #5
                      VRPipe::Step->get(name => 'vcf_index'),            #6
                      VRPipe::Step->get(name => 'vcf_concat'),           #7
                      VRPipe::Step->get(name => 'vcf_index'),            #8
                     ],

                     [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3, from_key => 'chunked_regions_file', to_key => 'chunked_regions_file'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 4, from_key => 'chunked_regions_file', to_key => 'chunked_regions_file'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'bcf_files', to_key => 'bcf_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 5, from_key => 'bcf_files', to_key => 'bcf_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'sample_ploidy_files', to_key => 'samples_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'vcf_files', to_key => 'vcf_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 7, from_key => 'vcf_files', to_key => 'vcf_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 8, from_key => 'concat_vcf', to_key => 'vcf_files'),
                     ],

                     [ VRPipe::StepBehaviourDefiner->new(after_step => 4, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1),
                       VRPipe::StepBehaviourDefiner->new(after_step => 6, behaviour => 'delete_outputs', act_on_steps => [4], regulated_by => 'cleanup', default_regulation => 1),
                       VRPipe::StepBehaviourDefiner->new(after_step => 6, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'remove_bcf', default_regulation => 0),
                       VRPipe::StepBehaviourDefiner->new(after_step => 7, behaviour => 'delete_outputs', act_on_steps => [5,6], regulated_by => 'cleanup', default_regulation => 1),
                     ]);
    }
}

1;