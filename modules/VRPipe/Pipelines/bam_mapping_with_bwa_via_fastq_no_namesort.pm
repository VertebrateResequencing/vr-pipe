
=head1 NAME

VRPipe::Pipelines::bam_mapping_with_bwa_via_fastq_no_namesort - a pipeline

=head1 DESCRIPTION

Maps reads in a bam file datasource to a reference genome using bwa fastq
alignment. For this the bams are first converted to fastq format, then
converted back after the alignment.

Since the bam_to_fastq step now uses the bam2fastq exe which does not require
name sorted bams, this pipeline avoids the expensive bam_name_sort step.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Pipelines::bam_mapping_with_bwa_via_fastq_no_namesort with VRPipe::PipelineRole {
    method name {
        return 'bam_mapping_with_bwa_via_fastq_no_namesort';
    }
    
    method description {
        return 'Map reads in bam files to a reference genome with bwa fastq alignment';
    }
    
    method step_names {
        (
            'sequence_dictionary',   #1
            'bwa_index',             #2
            'bam_metadata',          #3
            'bam_to_fastq',          #4
            'fastq_split',           #5
            'bwa_aln_fastq',         #6
            'bwa_sam',               #7
            'sam_to_fixed_bam',      #8
            'bam_merge_lane_splits', #9
            'bam_index',             #10
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 3,  to_key   => 'bam_files' },
            { from_step => 0, to_step => 4,  to_key   => 'bam_files' },
            { from_step => 4, to_step => 5,  from_key => 'fastq_files', to_key => 'fastq_files' },
            { from_step => 5, to_step => 6,  from_key => 'split_fastq_files', to_key => 'fastq_files' },
            { from_step => 6, to_step => 7,  from_key => 'bwa_sai_files', to_key => 'sai_files' },
            { from_step => 5, to_step => 7,  from_key => 'split_fastq_files', to_key => 'fastq_files' },
            { from_step => 7, to_step => 8,  from_key => 'bwa_sam_files', to_key => 'sam_files' },
            { from_step => 8, to_step => 9,  from_key => 'fixed_bam_files', to_key => 'bam_files' },
            { from_step => 1, to_step => 9,  from_key => 'reference_dict', to_key => 'dict_file' },
            { from_step => 9, to_step => 10, from_key => 'merged_lane_bams', to_key => 'bam_files' }
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 7, behaviour => 'delete_outputs', act_on_steps => [4, 5, 6], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 10, behaviour => 'delete_outputs', act_on_steps => [7, 8], regulated_by => 'cleanup', default_regulation => 1 }
        );
    }
}

1;
