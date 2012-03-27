use VRPipe::Base;

class VRPipe::Pipelines::bam_qc_graphs_and_stats with VRPipe::PipelineRole {
    method name {
        return 'bam_qc_graphs_and_stats';
    }
    method _num_steps {
        return 3;
    }
    method description {
        return 'Generates graphs and stats to describe bam files for quality-checking purposes.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'fasta_gc_stats'),#1
                  VRPipe::Step->get(name => 'bamcheck'),#2
                  VRPipe::Step->get(name => 'plot_bamcheck'),#3
                  VRPipe::Step->get(name => 'bamcheck_stats_output'),#4
                  ],

                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3, from_key => 'fasta_gc_stats_file', to_key => 'fasta_gc_stats_file'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'bamcheck_files', to_key => 'bamcheck_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 4, to_key => 'bam_files'),
                   ],

                 [ ]
                 );
    }
}

1;