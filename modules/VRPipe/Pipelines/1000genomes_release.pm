use VRPipe::Base;

class VRPipe::Pipelines::1000genomes_release with VRPipe::PipelineRole {
    method name {
        return '1000genomes_release';
    }
    method _num_steps {
        return 6;
    }
    method description {
        return 'Create 1000 genomes release files';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
            return ([ VRPipe::Step->get(name => 'dcc_metadata'),#1
                      VRPipe::Step->get(name => 'bam_index'),#2
                      VRPipe::Step->get(name => 'bam_stats'),#3
                      VRPipe::Step->get(name => 'md5_file_production'),#4
                      VRPipe::Step->get(name => 'md5_file_production'),#5
                      VRPipe::Step->get(name => 'md5_file_production'),#6
                      ],

                     [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'dcc_ready_bam_files', to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3, from_key => 'dcc_ready_bam_files', to_key => 'bam_files'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 4, from_key => 'dcc_ready_bam_files', to_key => 'md5_file_input'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 5, from_key => 'bai_files', to_key => 'md5_file_input'),
                       VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 6, from_key => 'bas_files', to_key => 'md5_file_input'),
                      ],

                     [ 
                     ]);
    }
}

1;