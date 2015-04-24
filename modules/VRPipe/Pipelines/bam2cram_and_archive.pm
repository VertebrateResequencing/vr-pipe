use VRPipe::Base;

class VRPipe::Pipelines::bam2cram_and_archive with VRPipe::PipelineRole {
    method name {
        return 'bam2cram_and_archive';
    }
    
    method description {
        return 'Convert BAM to CRAM and archive';
    }
    
    method step_names {
        (
            'bam_to_cram',   #1
            'archive_files', #2
            'cram_index',    #3
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 1, to_step => 1, from_key => 'cram_files', to_key => 'file' },
            { from_step => 1, to_step => 3, from_key => 'cram_files', to_key => 'cram_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 1, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_original_bams', default_regulation => 0 },
        );
    }
}

1;
