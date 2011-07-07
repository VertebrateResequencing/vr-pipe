use VRPipe::Base;

class VRPipe::Pipelines::mapping with VRPipe::PipelineRole {
    method name {
        return 'mapping';
    }
    method _num_steps {
        return 5;
    }
    method description {
        return 'Map reads in fastq files to a reference genome';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'fastq_metadata'),
                 #VRPipe::Step->get(name => 'fastq_split'),
                 #VRPipe::Step->get(name => 'fastq_map'),
                 #VRPipe::Step->get(name => 'bam_merge'),
                 #VRPipe::Step->get(name => 'bam_stats')
                 ],
                 [ {before_step_number => 1, adaptor_hash => { fastq_files => 'data_element' }},
                   {before_step_number => 2, adaptor_hash => { fastq_files => 'fastq_files_with_metadata' }},
                   {before_step_number => 3, adaptor_hash => { fastq_files => 'split_fastq_files' }},
                   {before_step_number => 4, adaptor_hash => { bam_files => 'mapped_bam_files' }},
                   {before_step_number => 5, adaptor_hash => { bam_files => 'merged_bam_files' }}]);
    }
}

1;