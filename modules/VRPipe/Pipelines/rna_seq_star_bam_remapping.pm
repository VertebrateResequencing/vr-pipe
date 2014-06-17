
=head1 NAME

VRPipe::Pipelines::rna_seq_star_bam_remapping - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Pipelines::rna_seq_star_bam_remapping with VRPipe::PipelineRole {
    method name {
        return 'rna_seq_star_bam_remapping';
    }
    
    method description {
        return 'Realign RNA-seq bam files with fast STAR algorithm.';
    }
    
    method step_names {
        (
            'sequence_dictionary', #1
            'star_buildgenome',    #2
            'bam_metadata',        #3
            'bam_to_fastq',        #4
            'star_map_fastq',      #5
            'sam_to_fixed_bam',    #6
            'bam_reheader',        #7
            'bam_index',           #8
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 3, to_key   => 'bam_files' },
            { from_step => 0, to_step => 4, to_key   => 'bam_files' },
            { from_step => 4, to_step => 5, from_key => 'fastq_files', to_key => 'fastq_files' },
            { from_step => 5, to_step => 6, from_key => 'star_sam_files', to_key => 'sam_files' },
            { from_step => 1, to_step => 7, from_key => 'reference_dict', to_key => 'dict_file' },
            { from_step => 6, to_step => 7, from_key => 'fixed_bam_files', to_key => 'bam_files' },
            { from_step => 7, to_step => 8, from_key => 'headed_bam_files', to_key => 'bam_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 8, behaviour => 'delete_outputs', act_on_steps => [4, 5, 6], regulated_by => 'cleanup', default_regulation => 1 },
        );
    }
}

1;
