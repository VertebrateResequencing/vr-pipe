
=head1 NAME

VRPipe::Pipelines::fermikit_calling - a pipeline

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

class VRPipe::Pipelines::fermikit_calling with VRPipe::PipelineRole {
    method name {
        return 'fermikit_calling';
    }
    
    method description {
        return 'Run the fermikit unitig calling pipeline.';
    }
    
    method step_names {
        (
            'fasta_index',           #1
            'sequence_dictionary',   #2
            'bwa_index',             #3
            'bwa_mem_align_unitigs', #4
            'bam_sort',              #5
            'unitig_pileup_calling', #6
            'unitig_abreak_calling', #7
            # 'filter_unitig_pileup_calls', #8
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 4, to_key   => 'fastq_files' },
            { from_step => 2, to_step => 4, from_key => 'reference_dict', to_key => 'dict_file' },
            { from_step => 4, to_step => 5, from_key => 'bwa_mem_bam_files', to_key => 'bam_files' },
            { from_step => 5, to_step => 6, from_key => 'coord_sorted_bam_files', to_key => 'aln_files' },
            { from_step => 4, to_step => 7, from_key => 'bwa_mem_bam_files', to_key => 'name_sorted_bam_files' },
            # { from_step => 6, to_step => 8, from_key => 'unitig_pileup_vcf_files', to_key => 'vcf_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 4, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'delete_input_files',     default_regulation => 0 },
            { after_step => 6, behaviour => 'delete_outputs', act_on_steps => [5], regulated_by => 'remove_alignment_files', default_regulation => 0 },
            { after_step => 7, behaviour => 'delete_outputs', act_on_steps => [4], regulated_by => 'remove_alignment_files', default_regulation => 0 },
            { after_step => 8, behaviour => 'delete_outputs', act_on_steps => [6], regulated_by => 'remove_unfiltered_vcfs', default_regulation => 0 },
        );
    }
}

1;
