use VRPipe::Base;

class VRPipe::Pipelines::graphs_and_stats_wgs with VRPipe::PipelineRole {
    method name {
        return 'graphs_and_stats_wgs';
    }
    method _num_steps {
        return 3;
        #return 2;
    }
    method description {
        return 'Performs stats and graph analysis for whole genome study sample bam files.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'qc_stats_bamcheck_wgs'),#1
                  VRPipe::Step->get(name => 'qc_stats_bamcheck_rmdup_wgs'),#2
                  VRPipe::Step->get(name => 'qc_plots_wgs'),#3
                  ],

                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'bam_files_with_metadata', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'bamcheck_files', to_key => 'bamcheck_files'),
                   ],

                 [ 
                  VRPipe::StepBehaviourDefiner->new(after_step => 3, behaviour => 'delete_outputs', act_on_steps => [2], regulated_by => 'cleanup', default_regulation => 0),
                 ]
                 );
    }
}

1;