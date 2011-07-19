use VRPipe::Base;

class VRPipe::Pipelines::fastq_mapping_with_bwa with VRPipe::PipelineRole {
    method name {
        return 'fastq_mapping_with_bwa';
    }
    method _num_steps {
        return 8;
    }
    method description {
        return 'Map reads in fastq files to a reference genome with bwa';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'fastq_metadata'),
                  VRPipe::Step->get(name => 'fastq_split'),
                  VRPipe::Step->get(name => 'bwa_index'),
                  VRPipe::Step->get(name => 'bwa_aln_fastq'),
                  VRPipe::Step->get(name => 'bwa_sam'),
                 #VRPipe::Step->get(name => 'sam_to_fixed_bam'),
                 #VRPipe::Step->get(name => 'bam_merge_lane_splits'),
                 #VRPipe::Step->get(name => 'bam_stats')
                 ],
                 [ {before_step_number => 1, adaptor_hash => { fastq_files => 'data_element' }},
                   {before_step_number => 2, adaptor_hash => { fastq_files => 'fastq_files_with_metadata' }},
                   {before_step_number => 4, adaptor_hash => { fastq_files => 'split_fastq_files' }},
                   {before_step_number => 5, adaptor_hash => { fastq_files => 'split_fastq_files', sai_files => 'bwa_sai_files' }},
                   {before_step_number => 6, adaptor_hash => { sam_files => 'bwa_sam_files' }},
                   {before_step_number => 7, adaptor_hash => { bam_files => 'fixed_bam_files' }},
                   {before_step_number => 8, adaptor_hash => { bam_files => 'merged_lane_bams' }}]);
    }
}

1;