use VRPipe::Base;

class VRPipe::Pipelines::graphs_and_stats_exome with VRPipe::PipelineRole {
    method name {
        return 'graphs_and_stats_exome';
    }
    method _num_steps {
        return 2;
    }
    method description {
        return 'Performs stats and graph analysis for exome study sample bam files.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'qc_stats_exome'),#1
                  VRPipe::Step->get(name => 'qc_plots_exome'),#2
                  ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'stats_dump_files', to_key => 'stats_dump_files'),
                   ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 0),
                 ]
                 );
    }
}

1;