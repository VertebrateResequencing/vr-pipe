use VRPipe::Base;

class VRPipe::Pipelines::bam_import_from_irods_and_vrtrack_qc_wgs with VRPipe::PipelineRole {
    method name {
        return 'bam_import_from_irods_and_vrtrack_qc_wgs';
    }
    method _num_steps {
        return 6;
    }
    method description {
        return 'Copy bam files stored in iRODs to local disc, generate QC stats and graphs and update the VRTrack db. For whole genome shotgun sequence only.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'irods_get_files_by_basename'),#1
                  #VRPipe::Step->get(name => 'bamcheck'),#2
                  #VRPipe::Step->get(name => 'plot_bamcheck'),#3
                  #VRPipe::Step->get(name => 'bamcheck_rmdup'),#4
                  #VRPipe::Step->get(name => 'bamcheck_stats_output'),#5
                  #VRPipe::Step->get(name => 'vrtrack_update_mapstats'),#6
                  ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'basenames'),
                   #VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'local_files', to_key => 'bam_files'),
                   #VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'bamcheck_files', to_key => 'bamcheck_files'),
                   #VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 4, from_key => 'local_files', to_key => 'bam_files'),
                   #VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'bam_files_with_metadata', to_key => 'bam_files'),
                   #VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 6, from_key => 'bam_files_with_metadata', to_key => 'bam_files')
                  ],
                 
                 [ ]);
    }
}

1;