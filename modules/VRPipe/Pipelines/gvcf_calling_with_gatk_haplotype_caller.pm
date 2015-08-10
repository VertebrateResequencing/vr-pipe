
=head1 NAME

VRPipe::Pipelines::gvcf_calling_with_gatk_haplotype_caller - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Pipelines::gvcf_calling_with_gatk_haplotype_caller with VRPipe::PipelineRole {
    method name {
        return 'gvcf_calling_with_gatk_haplotype_caller';
    }
    
    method description {
        return 'Run GATK to genotype gVCF files. Use this pipeline for up to 200 samples.';
    }
    
    method step_names {
        (
            'gatk_genotype_gvcfs_with_genome_chunking', #1
            'bcftools_concat',                          #2
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'gvcf_files' },
            { from_step => 0, to_step => 1, to_key   => 'gvcf_index_files' },
            { from_step => 1, to_step => 2, from_key => 'genotype_gvcf_file', to_key => 'vcf_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 1, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'delete_input_gvcfs', default_regulation => 0 },
            { after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'delete_gvcf_chunks', default_regulation => 1 },
        );
    }

}

1;
