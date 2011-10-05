use VRPipe::Base;

class VRPipe::Pipelines::improvement_test_pipeline with VRPipe::PipelineRole {
    method name {
        return 'improvement_test_pipeline';
    }
    method _num_steps {
        return 10;
    }
    method description {
        return 'a pipeline for testing bam improvement';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
            return ([ VRPipe::Step->get(name => 'test_import_bams'), #1
                      VRPipe::Step->get(name => 'sequence_dictionary'),#2
                      VRPipe::Step->get(name => 'bam_index'),#3
                      VRPipe::Step->get(name => 'gatk_target_interval_creator'),#4
                      VRPipe::Step->get(name => 'bam_realignment_around_known_indels'),#5
                      VRPipe::Step->get(name => 'bam_fix_mates'),#6 # should be able to skip this -- need to check
                      VRPipe::Step->get(name => 'bam_index'),#7
                      VRPipe::Step->get(name => 'bam_count_covariates'),#8
                      VRPipe::Step->get(name => 'bam_recalibrate_quality_scores'),#9
                      VRPipe::Step->get(name => 'bam_calculate_bq'),#10
                      ],

                     [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3, from_key => 'local_bam_files', to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 5, from_key => 'local_bam_files', to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'intervals_file', to_key => 'intervals_file'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'realigned_bam_files', to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'fixmate_bam_files', to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 8, from_key => 'fixmate_bam_files', to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 8, to_step => 9, from_key => 'bam_recalibration_files', to_key => 'bam_recalibration_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 9, from_key => 'fixmate_bam_files', to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 9, to_step => 10, from_key => 'recalibrated_bam_files', to_key => 'bam_files'),
                      ],

                     [ VRPipe::StepBehaviourDefiner->new(after_step => 6, behaviour => 'delete_outputs', act_on_steps => [5], regulated_by => 'cleanup', default_regulation => 1),
                       VRPipe::StepBehaviourDefiner->new(after_step => 9, behaviour => 'delete_outputs', act_on_steps => [6,7], regulated_by => 'cleanup', default_regulation => 1),
                       VRPipe::StepBehaviourDefiner->new(after_step => 10, behaviour => 'delete_outputs', act_on_steps => [8,9], regulated_by => 'cleanup', default_regulation => 1),
                      ]);
    }
}

1;