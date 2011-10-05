use VRPipe::Base;

class VRPipe::Pipelines::1000genomes_454_mapping_with_improvement with VRPipe::PipelineRole {
    method name {
        return '1000genomes_454_mapping_with_improvement';
    }
    method _num_steps {
        return 17;
    }
    method description {
        return 'Map (with improvement) LS454 reads in fastq files on the DCC ftp site to a reference genome with smalt';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'sequence_dictionary'),#1
                  VRPipe::Step->get(name => 'smalt_index'),#2
                  VRPipe::Step->get(name => 'fastq_import'),#3
                  VRPipe::Step->get(name => 'fastq_metadata'),#4
                  VRPipe::Step->get(name => 'fastq_split'),#5
                  VRPipe::Step->get(name => 'smalt_map_to_sam'),#6
                  VRPipe::Step->get(name => 'sam_to_fixed_bam'),#7
                  VRPipe::Step->get(name => 'bam_merge_lane_splits'),#8
                  VRPipe::Step->get(name => 'bam_index'),#9
                  VRPipe::Step->get(name => 'gatk_target_interval_creator'),#10
                  VRPipe::Step->get(name => 'bam_realignment_around_known_indels'),#11
                  VRPipe::Step->get(name => 'bam_index'),#12
                  VRPipe::Step->get(name => 'bam_count_covariates'),#13
                  VRPipe::Step->get(name => 'bam_recalibrate_quality_scores'),#14
                  VRPipe::Step->get(name => 'bam_calculate_bq'),#15
                  VRPipe::Step->get(name => 'bam_reheader'),#16
                  VRPipe::Step->get(name => 'bam_stats')#17
                  ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'local_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'fastq_files_with_metadata', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'split_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 6, from_key => 'smalt_index_binary_files', to_key => 'index_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'smalt_sam_files', to_key => 'sam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 8, from_key => 'fixed_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 8, from_key => 'reference_dict', to_key => 'dict_file'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 8, to_step => 9, from_key => 'merged_lane_bams', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 10, to_step => 11, from_key => 'intervals_file', to_key => 'intervals_file'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 8, to_step => 11, from_key => 'merged_lane_bams', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 11, to_step => 12, from_key => 'realigned_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 11, to_step => 13, from_key => 'realigned_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 13, to_step => 14, from_key => 'bam_recalibration_files', to_key => 'bam_recalibration_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 11, to_step => 14, from_key => 'realigned_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 14, to_step => 15, from_key => 'recalibrated_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 15, to_step => 16, from_key => 'bq_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 16, from_key => 'reference_dict', to_key => 'dict_file'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 16, to_step => 17, from_key => 'headed_bam_files', to_key => 'bam_files') 
                  ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 6, behaviour => 'delete_outputs', act_on_steps => [3, 5], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 7, behaviour => 'delete_outputs', act_on_steps => [6], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 8, behaviour => 'delete_outputs', act_on_steps => [7], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 11, behaviour => 'delete_outputs', act_on_steps => [8,9], regulated_by => 'remove_mapped_bams', default_regulation => 0),
                   VRPipe::StepBehaviourDefiner->new(after_step => 14, behaviour => 'delete_outputs', act_on_steps => [11,12,13], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 15, behaviour => 'delete_outputs', act_on_steps => [14], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 16, behaviour => 'delete_outputs', act_on_steps => [15], regulated_by => 'cleanup', default_regulation => 1) 
                  ]);
    }
}

1;