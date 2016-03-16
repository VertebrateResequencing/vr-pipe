
=head1 NAME

VRPipe::Pipelines::bam_calling_with_freebayes_with_vep - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2016 Genome Research Limited.

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

class VRPipe::Pipelines::bam_calling_with_freebayes_with_vep with VRPipe::PipelineRole {
    method name {
        return 'bam_calling_with_freebayes_with_vep';
    }
    
    method description {
        return 'Run freebayes split over the genome and concatenate the results. Also run Variant Effect Predictor (VEP).';
    }
    
    method step_names {
        (
            'freebayes_with_genome_chunking', #1
            'vep_annotate',                   #2
            'bcf_concat',                     #3
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 0, to_step => 1, to_key   => 'bam_index_files' },
            { from_step => 1, to_step => 2, from_key => 'freebayes_variant_files', to_key => 'vcf_files' },
            { from_step => 2, to_step => 3, from_key => 'vep_annot_vcf_files', to_key => 'vcf_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 1, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0 },
            { after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup',           default_regulation => 1 },
            { after_step => 3, behaviour => 'delete_outputs', act_on_steps => [2], regulated_by => 'cleanup',           default_regulation => 1 },
        );
    }
}

1;
