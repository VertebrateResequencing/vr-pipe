=head1 NAME

VRPipe::Pipelines::snp_calling_chunked_mpileup_vcf - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::Pipelines::snp_calling_chunked_mpileup_vcf with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_chunked_mpileup_vcf';
    }
    method _num_steps {
        return 3;
    }
    method description {
        return 'Run samtools mpileup with bam chunking to generate vcf files';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
					VRPipe::Step->get(name => 'chunk_genomic_region'),
					VRPipe::Step->get(name => 'mpileup_vcf'),
					VRPipe::Step->get(name => 'vcf_concat'),
				],
   	            [ 
					VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'),
					VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'chunked_regions_file', to_key => 'chunked_regions_file'), 
					VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'vcf_files', to_key => 'vcf_files'), 
				],
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 3, behaviour => 'delete_outputs', act_on_steps => [2], regulated_by => 'cleanup', default_regulation => 1) ]	
				);
    }
}

1;
