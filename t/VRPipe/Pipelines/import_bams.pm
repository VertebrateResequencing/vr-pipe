use VRPipe::Base;

class VRPipe::Pipelines::import_bams with VRPipe::PipelineRole {
    method name {
        return 'import_bams';
    }
    
    method description {
        return 'a pipeline for importing bams, getting their metadata';
    }
    
    method step_names {
        (
            'test_import_bams', #1
            'bam_metadata',     #2
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 1, to_step => 2, from_key => 'local_bam_files', to_key => 'bam_files' }
        );
    }
}

1;
