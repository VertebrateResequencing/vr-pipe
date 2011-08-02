use VRPipe::Base;

class VRPipe::Pipelines::1000genomes_illumina_mapping_with_improvement with VRPipe::PipelineRole {
    method name {
        return '1000genomes_illumina_mapping_with_improvement';
    }
    method _num_steps {
        return 15;
    }
    method description {
        return 'Map (with improvement) Illumina reads in fastq files on the DCC ftp site to a reference genome with bwa';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'sequence_dictionary'),
                  VRPipe::Step->get(name => 'bwa_index'),
                  VRPipe::Step->get(name => 'fastq_import'),
                  VRPipe::Step->get(name => 'fastq_metadata'),
                  VRPipe::Step->get(name => 'fastq_split'),
                  VRPipe::Step->get(name => 'bwa_aln_fastq'),
                  VRPipe::Step->get(name => 'bwa_sam'),
                  VRPipe::Step->get(name => 'sam_to_fixed_bam'),
                  VRPipe::Step->get(name => 'bam_merge_lane_splits'),
                  
                  VRPipe::Step->get(name => 'bam_realignment_around_known_indels'),
                  VRPipe::Step->get(name => 'bam_fix_mates'),
                  VRPipe::Step->get(name => 'bam_recalibrate_quality_scores'),
                  VRPipe::Step->get(name => 'bam_calculate_bq'),
                  VRPipe::Step->get(name => 'bam_reheader'),
                  
                  VRPipe::Step->get(name => 'bam_stats') ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'local_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'fastq_files_with_metadata', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'split_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 7, from_key => 'split_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'bwa_sai_files', to_key => 'sai_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 8, from_key => 'bwa_sam_files', to_key => 'sam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 8, to_step => 9, from_key => 'fixed_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 9, from_key => 'reference_dict', to_key => 'dict_file'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 9, to_step => 10, from_key => 'merged_lane_bams', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 10, to_step => 11, from_key => 'realigned_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 11, to_step => 12, from_key => 'fixmate_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 12, to_step => 13, from_key => 'recalibrated_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 13, to_step => 14, from_key => 'bq_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 14, from_key => 'reference_dict', to_key => 'dict_file'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 14, to_step => 15, from_key => 'headed_bam_files', to_key => 'bam_files') ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 7, behaviour => 'delete_outputs', act_on_steps => [3, 5, 6], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 8, behaviour => 'delete_outputs', act_on_steps => [7], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 9, behaviour => 'delete_outputs', act_on_steps => [8], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 10, behaviour => 'delete_outputs', act_on_steps => [9], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 11, behaviour => 'delete_outputs', act_on_steps => [10], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 12, behaviour => 'delete_outputs', act_on_steps => [11], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 13, behaviour => 'delete_outputs', act_on_steps => [12], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;