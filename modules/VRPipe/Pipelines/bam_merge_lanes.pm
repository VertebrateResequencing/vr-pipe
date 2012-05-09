use VRPipe::Base;

class VRPipe::Pipelines::bam_merge_lanes with VRPipe::PipelineRole {
    method name {
        return 'bam_merge_lanes';
    }
    method _num_steps {
        return 3;
    }
    method description {
        return 'Tag strip, merge and mark duplicates';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
            return ([ VRPipe::Step->get(name => 'bam_strip_tags'),#1
                      VRPipe::Step->get(name => 'bam_merge'),#2
                      VRPipe::Step->get(name => 'bam_mark_duplicates'),#3
                      ],

                     [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'tag_stripped_bam_files', to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'merged_bam_files', to_key => 'bam_files'),
                      ],

                     [ VRPipe::StepBehaviourDefiner->new(after_step => 1, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0),
                       VRPipe::StepBehaviourDefiner->new(after_step => 3, behaviour => 'delete_outputs', act_on_steps => [1,2], regulated_by => 'cleanup', default_regulation => 1)
                     ]);
    }
}

1;