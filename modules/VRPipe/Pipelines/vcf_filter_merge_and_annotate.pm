use VRPipe::Base;

class VRPipe::Pipelines::vcf_filter_merge_and_annotate with VRPipe::PipelineRole {
    method name {
        return 'vcf_filter_merge_and_annotate';
    }
    method _num_steps {
        return 4;
    }
    method description {
        return 'Filter then merge VCF files and then annotate them, including consequences';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'vcf_filter'),
                  VRPipe::Step->get(name => 'vcf_merge'),
                  VRPipe::Step->get(name => 'vcf_annotate'),
                  VRPipe::Step->get(name => 'vcf_consequences') ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'vcf_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'filtered_vcf', to_key => 'vcf_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'merged_vcf', to_key => 'vcf_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'annotated_vcf', to_key => 'vcf_files') ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 3, behaviour => 'delete_outputs', act_on_steps => [2], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 4, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;
