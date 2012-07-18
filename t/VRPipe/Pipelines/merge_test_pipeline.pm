use VRPipe::Base;

class VRPipe::Pipelines::merge_test_pipeline with VRPipe::PipelineRole {
    method name {
        return 'merge_test_pipeline';
    }
    
    method _num_steps {
        return 5;
    }
    
    method description {
        return 'a pipeline for testing bam merging';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
             VRPipe::Step->get(name => 'test_import_bams'),    #1
             VRPipe::Step->get(name => 'bam_metadata'),        #2
             VRPipe::Step->get(name => 'bam_strip_tags'),      #3
             VRPipe::Step->get(name => 'bam_merge'),           #4
             VRPipe::Step->get(name => 'bam_mark_duplicates'), #5
            ],
            
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'local_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3, from_key => 'local_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'tag_stripped_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'merged_bam_files', to_key => 'bam_files'),],
            
            [VRPipe::StepBehaviourDefiner->new(after_step => 3, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1), VRPipe::StepBehaviourDefiner->new(after_step => 5, behaviour => 'delete_outputs', act_on_steps => [3, 4], regulated_by => 'cleanup', default_regulation => 1)]);
    }
}

1;
