
=head1 NAME

VRPipe::Pipelines::sample_bam_remapping_with_bwa_mem_nosplit - a pipeline

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

class VRPipe::Pipelines::sample_bam_remapping_with_bwa_mem_nosplit with VRPipe::PipelineRole {
    method name {
        return 'sample_bam_remapping_with_bwa_mem_nosplit';
    }
    
    method description {
        return 'Remap BAM or CRAM files for a sample using the bwa mem algorithm. Each sample should have one BAM or CRAM file per-readgroup.';
    }
    
    method step_names {
        (
            'fasta_index',                                      #1
            'sequence_dictionary',                              #2
            'bwa_index',                                        #3
            'bamtofastq',                                       #4 filter out SECONDARY,SUPPLEMENTARY,QCFAIL reads, rescue orphans 'gz=1 exclude=SECONDARY,SUPPLEMENTARY,QCFAIL level=1'
            'bwa_mem_to_bam',                                   #5 this will add sequence dictionary and @RG line, biobambam/bamsort inputformat=sam fixmates=1 adddupmarksupport=1 level=1
            'biobambam_bammerge_and_streaming_mark_duplicates', #6
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 4, to_key   => 'bam_files' },
            { from_step => 4, to_step => 5, from_key => 'fastq_files', to_key => 'fastq_files' },
            { from_step => 2, to_step => 5, from_key => 'reference_dict', to_key => 'dict_file' },
            { from_step => 5, to_step => 6, from_key => 'bam_files', to_key => 'bam_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 4, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0 },
            { after_step => 5, behaviour => 'delete_outputs', act_on_steps => [4], regulated_by => 'cleanup',           default_regulation => 1 },
            { after_step => 6, behaviour => 'delete_outputs', act_on_steps => [5], regulated_by => 'cleanup',           default_regulation => 1 },
        );
    }
}

1;
