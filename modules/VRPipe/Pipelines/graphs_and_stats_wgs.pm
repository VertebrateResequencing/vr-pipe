use VRPipe::Base;

class VRPipe::Pipelines::graphs_and_stats_wgs with VRPipe::PipelineRole {
    method name {
        return 'graphs_and_stats_wgs';
    }
    method _num_steps {
        return 3;
    }
    method description {
        return 'Performs stats and graph analysis for whole genome study sample bam files.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'bamcheck'),#1
                  VRPipe::Step->get(name => 'plot_bamcheck'),#2
                  VRPipe::Step->get(name => 'bamcheck_stats_output'),#3
                  ],

                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'bamcheck_files', to_key => 'bamcheck_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'bam_files'),
                   ],

                 [ ]
                 );
    }
}

1;