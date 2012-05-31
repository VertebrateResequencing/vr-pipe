use VRPipe::Base;

class VRPipe::Pipelines::retroseq_analysis with VRPipe::PipelineRole {
    method name {
        return 'retroseq_analysis';
    }
    method _num_steps {
        return 2;
    }
    method description {
        return 'Run retroseq genotyping of transposable elements from short read alignments';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ 
                  VRPipe::Step->get(name => 'retroseq_discover'), 
                  VRPipe::Step->get(name => 'retroseq_call'),
                ],
                 
                [ 
                  VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                  VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'),
                  VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'rseq_bed', to_key => 'rseq_bed'),
                ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;
