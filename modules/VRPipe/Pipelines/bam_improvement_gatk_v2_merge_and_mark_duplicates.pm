
=head1 NAME

VRPipe::Pipelines::bam_improvement_gatk_v2_merge_and_mark_duplicates - a
pipeline

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

class VRPipe::Pipelines::bam_improvement_gatk_v2_merge_and_mark_duplicates with VRPipe::PipelineRole {
    method name {
        return 'bam_improvement_gatk_v2_merge_and_mark_duplicates';
    }
    
    method description {
        return 'Runs the new versions of GATK (>v2.x) on input bam files to realign the reads around known indels and recalibrate the base quality scores. The input bam files must be grouped in the datasource (e.g. by library). After bam improvement using GATK, the resulting alignment files are merged together and duplicated reads are marked using biobambam to produce one bam per element. Do not set --disable_bam_indexing in any of the GATK step options';
    }
    
    method step_names {
        (
            'samtools_add_readgroup',                       #1
            'gatk_target_interval_creator',                 #2
            'bam_realignment_around_known_indels',          #3
            'gatk_base_recalibrator',                       #4
            'gatk_print_reads_with_bqsr',                   #5
            'samtools_merge_and_streaming_mark_duplicates', #6
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 1, to_step => 3, from_key => 'rg_added_bam_files', to_key => 'bam_files' },
            { from_step => 1, to_step => 3, from_key => 'rg_added_index_files', to_key => 'bai_files' },
            { from_step => 2, to_step => 3, from_key => 'intervals_file', to_key => 'intervals_file' },
            { from_step => 3, to_step => 4, from_key => 'realigned_bam_files', to_key => 'bam_files' },
            { from_step => 3, to_step => 4, from_key => 'realigned_bam_index_files', to_key => 'bai_files' },
            { from_step => 3, to_step => 5, from_key => 'realigned_bam_files', to_key => 'bam_files' },
            { from_step => 3, to_step => 5, from_key => 'realigned_bam_index_files', to_key => 'bai_files' },
            { from_step => 4, to_step => 5, from_key => 'bam_recalibration_files', to_key => 'bam_recalibration_files' },
            { from_step => 5, to_step => 6, from_key => 'recalibrated_bam_files', to_key => 'bam_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 3, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0 },
            { after_step => 5, behaviour => 'delete_outputs', act_on_steps => [1, 2, 3, 4], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 6, behaviour => 'delete_outputs', act_on_steps => [5], regulated_by => 'cleanup', default_regulation => 1 },
        );
    }
}

1;
