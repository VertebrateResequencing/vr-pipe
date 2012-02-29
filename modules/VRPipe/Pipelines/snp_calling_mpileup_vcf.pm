use VRPipe::Base;

class VRPipe::Pipelines::snp_calling_mpileup_vcf with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_mpileup_vcf';
    }
    method _num_steps {
        return 1;
    }
    method description {
        return 'Run samtools mpileup to generate vcf files';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'mpileup_vcf') ],
   	        [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files') ],
                [ ]);
    }
}

1;
