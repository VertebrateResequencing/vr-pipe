
=head1 NAME

VRPipe::Pipelines::sga_prepare_region_fastq - a pipeline

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

class VRPipe::Pipelines::sga_prepare_region_fastq with VRPipe::PipelineRole {
    method name {
        return 'sga_prepare_region_fastq';
    }
    
    method description {
        return 'Split BAMs by region, converts to fastq and preprocesses with sga.';
    }
    
    method step_names {
        (
            'bam_split_by_region', #1
            'bam_metadata',        #2
            'bam_to_fastq',        #3
            'sga_preprocess',      #4
            'fastq_metadata'       #5
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 1, to_step => 2, from_key => 'split_region_bam_files', to_key => 'bam_files' },
            { from_step => 1, to_step => 3, from_key => 'split_region_bam_files', to_key => 'bam_files' },
            { from_step => 3, to_step => 4, from_key => 'fastq_files', to_key => 'fastq_files' },
            { from_step => 4, to_step => 5, from_key => 'preprocessed_fastq_files', to_key => 'fastq_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 3, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 4, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'cleanup', default_regulation => 1 }
        );
    }
}

1;
