
=head1 NAME

VRPipe::Pipelines::sample_gvcf_with_gatk_haplotype_caller - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Pipelines::sample_gvcf_with_gatk_haplotype_caller with VRPipe::PipelineRole {
    method name {
        return 'sample_gvcf_with_gatk_haplotype_caller';
    }
    
    method description {
        return 'Run GATK Haplotype Caller for per-sample calling split over the genome; BAM and BAI files required as input';
    }
    
    method step_names {
        (
            'fasta_index',                                #1
            'sequence_dictionary',                        #2
            'gatk_haplotype_caller_with_genome_chunking', #3
            'bcftools_concat',                            #4
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 3, to_key   => 'bam_files' },
            { from_step => 0, to_step => 3, to_key   => 'bai_files' },
            { from_step => 0, to_step => 4, to_key   => 'sites_file' },
            { from_step => 3, to_step => 4, from_key => 'gatk_hc_vcf_files', to_key => 'vcf_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 3, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0 },
            { after_step => 4, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'cleanup',           default_regulation => 1 }
        );
    }
}

1;
