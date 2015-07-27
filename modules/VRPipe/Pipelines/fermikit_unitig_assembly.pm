
=head1 NAME

VRPipe::Pipelines::fermikit_unitig_assembly - a pipeline

=head1 DESCRIPTION

Pipeline implementing the unitig assembly part of the fermikit pipeline,
https://github.com/lh3/fermikit

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

class VRPipe::Pipelines::fermikit_unitig_assembly with VRPipe::PipelineRole {
    method name {
        return 'fermikit_unitig_assembly';
    }
    
    method description {
        return 'Run the fermikit unitig assembly.';
    }
    
    method step_names {
        (
            'bfc_error_correct',                #1
            'bfc_filter_error_corrected_reads', #2
            'fmd_index',                        #3
            'fermi2_assemble',                  #4
            'fermi2_simplify',                  #5
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 1, to_step => 2, from_key => 'bfc_error_corrected_reads', to_key => 'bfc_error_corrected_reads' },
            { from_step => 2, to_step => 3, from_key => 'bfc_filtered_error_corrected_reads', to_key => 'fastq_files' },
            { from_step => 3, to_step => 4, from_key => 'fmd_index_files', to_key => 'fmd_index_files' },
            { from_step => 4, to_step => 5, from_key => 'assembled_unitigs', to_key => 'assembled_unitigs' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 1, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0 },
            { after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup',           default_regulation => 1 },
            { after_step => 3, behaviour => 'delete_outputs', act_on_steps => [2], regulated_by => 'cleanup',           default_regulation => 1 },
            { after_step => 4, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'cleanup',           default_regulation => 1 },
            { after_step => 5, behaviour => 'delete_outputs', act_on_steps => [4], regulated_by => 'cleanup',           default_regulation => 1 },
        );
    }
}

1;
