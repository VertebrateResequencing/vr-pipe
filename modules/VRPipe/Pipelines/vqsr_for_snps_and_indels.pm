
=head1 NAME

VRPipe::Pipelines::vqsr_for_snps_and_indels - a pipeline

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

class VRPipe::Pipelines::vqsr_for_snps_and_indels with VRPipe::PipelineRole {
    method name {
        return 'vqsr_for_snps_and_indels';
    }
    
    method description {
        return 'Run VQSR on SNPs and indels.';
    }
    
    method step_names {
        (
            'vcf_index',                             #1
            'gatk_variant_recalibration_for_indels', #2
            'gatk_variant_recalibration_for_snps',   #3
            'gatk_apply_recalibration_for_indels',   #4
            'vcf_index',                             #5
            'gatk_apply_recalibration_for_snps',     #6
            'vcf_index',                             #7
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'vcf_files' },
            { from_step => 0, to_step => 2, to_key   => 'vcf_files' },
            { from_step => 0, to_step => 3, to_key   => 'vcf_files' },
            { from_step => 0, to_step => 4, to_key   => 'vcf_files' },
            { from_step => 2, to_step => 4, from_key => 'recalibration_file', to_key => 'recalibration_file' },
            { from_step => 2, to_step => 4, from_key => 'tranches_file', to_key => 'tranches_file' },
            { from_step => 4, to_step => 5, from_key => 'recalibrated_vcfs', to_key => 'vcf_files' },
            { from_step => 4, to_step => 6, from_key => 'recalibrated_vcfs', to_key => 'vcf_files' },
            { from_step => 3, to_step => 6, from_key => 'recalibration_file', to_key => 'recalibration_file' },
            { from_step => 3, to_step => 6, from_key => 'tranches_file', to_key => 'tranches_file' },
            { from_step => 6, to_step => 7, from_key => 'recalibrated_vcfs', to_key => 'vcf_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 5, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'remove_input_vcfs',          default_regulation => 0 },
            { after_step => 5, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'remove_input_vcfs',          default_regulation => 0 },
            { after_step => 8, behaviour => 'delete_outputs', act_on_steps => [2], regulated_by => 'remove_recalibration_files', default_regulation => 0 },
            { after_step => 6, behaviour => 'delete_outputs', act_on_steps => [4, 5], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 7, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'remove_recalibration_files', default_regulation => 0 },
        );
    }
}

1;
