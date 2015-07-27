
=head1 NAME

VRPipe::Pipelines::sample_fastq_mapping_with_bwa_mem - a pipeline

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

class VRPipe::Pipelines::sample_fastq_mapping_with_bwa_mem with VRPipe::PipelineRole {
    method name {
        return 'sample_fastq_mapping_with_bwa_mem';
    }
    
    method description {
        return 'Map reads in fastq files for a sample using the bwa mem algorithm';
    }
    
    method step_names {
        (
            'fasta_index',                                      #1
            'sequence_dictionary',                              #2
            'bwa_index',                                        #3
            'fastq_metadata',                                   #4
            'fastq_split',                                      #5
            'bwa_mem_to_bam',                                   #6 this will add sequence dictionary and @RG line, biobambam/bamsort inputformat=sam fixmates=1 adddupmarksupport=1 level=1
            'biobambam_bammerge_and_streaming_mark_duplicates', #7
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 4, to_key   => 'fastq_files' },
            { from_step => 0, to_step => 5, to_key   => 'fastq_files' },
            { from_step => 5, to_step => 6, from_key => 'split_fastq_files', to_key => 'fastq_files' },
            { from_step => 2, to_step => 6, from_key => 'reference_dict', to_key => 'dict_file' },
            { from_step => 6, to_step => 7, from_key => 'bam_files', to_key => 'bam_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 5, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'delete_input_fastq', default_regulation => 0 },
            { after_step => 6, behaviour => 'delete_outputs', act_on_steps => [5], regulated_by => 'cleanup',            default_regulation => 1 },
            { after_step => 7, behaviour => 'delete_outputs', act_on_steps => [6], regulated_by => 'cleanup',            default_regulation => 1 },
        );
    }
}

1;
