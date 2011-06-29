use VRPipe::Base;

class VRPipe::Pipelines::mapping with VRPipe::PipelineRole {
    method name {
        return 'mapping';
    }
    method _num_steps {
        return 4;
    }
    method description {
        return 'Map reads in fastq files to a reference genome';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return [ VRPipe::Step->get(name => 'fastq_split'),
                 VRPipe::Step->get(name => 'fastq_map'),
                 VRPipe::Step->get(name => 'bam_merge'),
                 VRPipe::Step->get(name => 'bam_stats') ];
    }
}

1;