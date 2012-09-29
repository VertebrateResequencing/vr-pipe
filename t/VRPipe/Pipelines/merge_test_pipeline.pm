use VRPipe::Base;

class VRPipe::Pipelines::merge_test_pipeline with VRPipe::PipelineRole {
    method name {
        return 'merge_test_pipeline';
    }
    
    method description {
        return 'a pipeline for testing bam merging';
    }
    
    method step_names {
        (
            'test_import_bams',    #1
            'bam_metadata',        #2
            'bam_strip_tags',      #3
            'bam_merge',           #4
            'bam_mark_duplicates', #5
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 1, to_step => 2, from_key => 'local_bam_files', to_key => 'bam_files' },
            { from_step => 1, to_step => 3, from_key => 'local_bam_files', to_key => 'bam_files' },
            { from_step => 3, to_step => 4, from_key => 'tag_stripped_bam_files', to_key => 'bam_files' },
            { from_step => 4, to_step => 5, from_key => 'merged_bam_files', to_key => 'bam_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 3, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 5, behaviour => 'delete_outputs', act_on_steps => [3, 4], regulated_by => 'cleanup', default_regulation => 1 }
        );
    }
}

1;
