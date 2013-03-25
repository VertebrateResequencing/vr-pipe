
=head1 NAME

VRPipe::Pipelines::bam_mapping_with_bwa - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Pipelines::bam_mapping_with_bwa with VRPipe::PipelineRole {
    method name {
        return 'bam_mapping_with_bwa';
    }
    
    method description {
        return 'Map reads in bam files to a reference genome with bwa. NB: this is unreliable and you may have better success with bam_mapping_with_bwa_via_fastq_no_namesort';
    }
    
    method step_names {
        (
            'sequence_dictionary', #1
            'bwa_index',           #2
            'bam_metadata',        #3
            'bam_name_sort',       #4
            'bwa_aln_bam',         #5
            'bwa_sam_using_bam',   #6
            'sam_to_fixed_bam',    #7
            'bam_reheader',        #8
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 3, to_key   => 'bam_files' },
            { from_step => 0, to_step => 4, to_key   => 'bam_files' },
            { from_step => 4, to_step => 5, from_key => 'name_sorted_bam_files', to_key => 'bam_files' },
            { from_step => 5, to_step => 6, from_key => 'bwa_sai_files', to_key => 'sai_files' },
            { from_step => 4, to_step => 6, from_key => 'name_sorted_bam_files', to_key => 'bam_files' },
            { from_step => 6, to_step => 7, from_key => 'bwa_sam_files', to_key => 'sam_files' },
            { from_step => 7, to_step => 8, from_key => 'fixed_bam_files', to_key => 'bam_files' },
            { from_step => 1, to_step => 8, from_key => 'reference_dict', to_key => 'dict_file' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 7, behaviour => 'delete_outputs', act_on_steps => [4, 5, 6], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 8, behaviour => 'delete_outputs', act_on_steps => [7], regulated_by => 'cleanup', default_regulation => 1 }
        );
    }
}

1;
