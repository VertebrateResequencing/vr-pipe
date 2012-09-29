
=head1 NAME

VRPipe::Pipelines::sga_merge_and_variant_calling - a pipeline

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

class VRPipe::Pipelines::sga_merge_and_variant_calling with VRPipe::PipelineRole {
    method name {
        return 'sga_merge_and_variant_calling';
    }
    
    method description {
        return 'Merge up to 16 fastq files and use sga to call variants: sga merge three times and sga graph-diff.';
    }
    
    method step_names {
        (
            'sga_permute_reference',             #1
            'sga_index_reference',               #2
            'sga_generate_sampled_suffix_array', #3
            'sga_merge',                         #4
            'sga_merge',                         #5
            'sga_merge',                         #6
            'sga_merge',                         #7
            'sga_reference_based_calling',       #8
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 4, to_key => 'fastq_files' },
            #            {from_step => 0, to_step => 4, to_key => 'bwt_files'},
            #            {from_step => 0, to_step => 4, to_key => 'sai_files'},
            { from_step => 1, to_step => 2, from_key => 'permuted_reference_fasta', to_key => 'reference_fasta' },
            { from_step => 1, to_step => 3, from_key => 'permuted_reference_fasta', to_key => 'reference_fasta' },
            { from_step => 4, to_step => 5, from_key => 'merged_fastq_files',       to_key => 'fastq_files' },
            #            {from_step => 4, to_step => 5, from_key => 'merged_bwt_files', to_key => 'bwt_files'},
            #            {from_step => 4, to_step => 5, from_key => 'merged_sai_files', to_key => 'sai_files'},
            { from_step => 5, to_step => 6, from_key => 'merged_fastq_files', to_key => 'fastq_files' },
            #            {from_step => 5, to_step => 6, from_key => 'merged_bwt_files', to_key => 'bwt_files'},
            #            {from_step => 5, to_step => 6, from_key => 'merged_sai_files', to_key => 'sai_files'},
            { from_step => 6, to_step => 7, from_key => 'merged_fastq_files', to_key => 'fastq_files' },
            #            {from_step => 6, to_step => 7, from_key => 'merged_bwt_files', to_key => 'bwt_files'},
            #            {from_step => 6, to_step => 7, from_key => 'merged_sai_files', to_key => 'sai_files'},
            { from_step => 1, to_step => 8, from_key => 'permuted_reference_fasta', to_key => 'reference_fasta' },
            { from_step => 7, to_step => 8, from_key => 'merged_fastq_files',       to_key => 'sga_indexed_variant_reads' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 4, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'remove_inputs',             default_regulation => 0 },
            { after_step => 5, behaviour => 'delete_outputs', act_on_steps => [4], regulated_by => 'cleanup',                   default_regulation => 1 },
            { after_step => 6, behaviour => 'delete_outputs', act_on_steps => [5], regulated_by => 'cleanup',                   default_regulation => 1 },
            { after_step => 7, behaviour => 'delete_outputs', act_on_steps => [6], regulated_by => 'cleanup',                   default_regulation => 1 },
            { after_step => 8, behaviour => 'delete_outputs', act_on_steps => [7], regulated_by => 'remove_final_merged_fastq', default_regulation => 0 }
        );
    }
}

1;
