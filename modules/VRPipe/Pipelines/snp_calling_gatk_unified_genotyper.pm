
=head1 NAME

VRPipe::Pipelines::snp_calling_gatk_unified_genotyper - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Pipelines::snp_calling_gatk_unified_genotyper with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_gatk_unified_genotyper';
    }
    
    method _num_steps {
        return 5;
    }
    
    method description {
        return 'Run gatk unified genotyper';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
                VRPipe::Step->get(name => 'fasta_index'),            #1 ## steps 1 and 2 to prevent runs on creating these files
                VRPipe::Step->get(name => 'sequence_dictionary'),    #2
                VRPipe::Step->get(name => 'bam_index'),              #3
                VRPipe::Step->get(name => 'gatk_unified_genotyper'), #4
                VRPipe::Step->get(name => 'vcf_index'),              #5
            ],
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 4, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'gatk_vcf_file', to_key => 'vcf_files'),],
            [],
        );
    }
}

1;
