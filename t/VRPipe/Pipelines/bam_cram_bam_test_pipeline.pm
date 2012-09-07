use VRPipe::Base;

class VRPipe::Pipelines::bam_cram_bam_test_pipeline with VRPipe::PipelineRole {
    method name {
        return 'bam_cram_bam_test_pipeline';
    }
    
    method _num_steps {
        return 7;
    }
    
    method description {
        return 'a pipeline for testing conversion of bam to cram and cram to bam';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
                VRPipe::Step->get(name => 'test_import_bams'), #1
                VRPipe::Step->get(name => 'fasta_index'),      #2
                VRPipe::Step->get(name => 'bam_metadata'),     #3
                VRPipe::Step->get(name => 'bam_index'),        #4
                VRPipe::Step->get(name => 'bam_to_cram'),      #5
                VRPipe::Step->get(name => 'cram_index'),       #6
                VRPipe::Step->get(name => 'cram_to_bam'),      #7
            ],
            
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3, from_key => 'local_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 4, from_key => 'local_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 5, from_key => 'local_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'bai_files', to_key => 'bai_files'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'cram_files', to_key => 'cram_files'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 7, from_key => 'cram_files', to_key => 'cram_files'), VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'cram_index_files', to_key => 'cram_index_files'),],
            
            []
        );
    }
}

1;
