
=head1 NAME

VRPipe::Pipelines::snp_calling_gatk_unified_genotyper_and_annotate - a pipeline

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

class VRPipe::Pipelines::snp_calling_gatk_unified_genotyper_and_annotate with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_gatk_unified_genotyper_and_annotate';
    }
    
    method description {
        return 'Run gatk unified genotyper followed by vcf-annotate';
    }
    
    method step_names {
        (
            'fasta_index',            #1 ## steps 1 and 2 to prevent runs on creating these files
            'sequence_dictionary',    #2
            'bam_index',              #3
            'gatk_unified_genotyper', #4
            'vcf_annotate',           #5
            'vcf_index',              #6
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 3, to_key   => 'bam_files' },
            { from_step => 0, to_step => 4, to_key   => 'bam_files' },
            { from_step => 0, to_step => 4, to_key   => 'sites_file' },
            { from_step => 4, to_step => 5, from_key => 'gatk_vcf_file', to_key => 'vcf_files' },
            { from_step => 5, to_step => 6, from_key => 'annotated_vcf', to_key => 'vcf_files' }
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 6, behaviour => 'delete_outputs', act_on_steps => [4], regulated_by => 'cleanup', default_regulation => 1 });
    }
}

1;
