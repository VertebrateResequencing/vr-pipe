use VRPipe::Base;

class VRPipe::Pipelines::bam_cram_bam_test_pipeline with VRPipe::PipelineRole {
    method name {
        return 'bam_cram_bam_test_pipeline';
    }
    
    method description {
        return 'a pipeline for testing conversion of bam to cram and cram to bam';
    }
    
    method step_names {
        (
            'test_import_bams', #1
            'fasta_index',      #2
            'bam_metadata',     #3
            'bam_index',        #4
            'bam_to_cram',      #5
            'cram_index',       #6
            'cram_to_bam',      #7
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 1, to_step => 3, from_key => 'local_bam_files', to_key => 'bam_files' },
            { from_step => 1, to_step => 4, from_key => 'local_bam_files', to_key => 'bam_files' },
            { from_step => 1, to_step => 5, from_key => 'local_bam_files', to_key => 'bam_files' },
            { from_step => 4, to_step => 5, from_key => 'bai_files', to_key => 'bai_files' },
            { from_step => 5, to_step => 6, from_key => 'cram_files', to_key => 'cram_files' },
            { from_step => 5, to_step => 7, from_key => 'cram_files', to_key => 'cram_files' },
            { from_step => 6, to_step => 7, from_key => 'cram_index_files', to_key => 'cram_index_files' },
        );
    }
}

1;
