
=head1 NAME

VRPipe::Pipelines::snp_calling_gatk_vcf - a pipeline

=head1 DESCRIPTION

Runs the GATK universal genotyper, followed by soft-filtering using the GATK VariantFiltration walker and recalibration of the variants using the GATK VariantRecalibrator walker, generating indexed VCF files from a BAM datasource.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Pipelines::snp_calling_gatk_vcf with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_gatk_vcf';
    }
    
    method _num_steps {
        return 7;
    }
    
    method description {
        return 'Run gatk universal genotyper and variant filtration/recalibration';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
             VRPipe::Step->get(name => 'gatk_genotype'),             #1
             VRPipe::Step->get(name => 'vcf_index'),                 #2
             VRPipe::Step->get(name => 'gatk_variant_filter'),       #3
             VRPipe::Step->get(name => 'vcf_index'),                 #4
             VRPipe::Step->get(name => 'gatk_recalibrate_variants'), #5
             VRPipe::Step->get(name => 'gatk_apply_recalibration'),  #6
             VRPipe::Step->get(name => 'vcf_index'),                 #7
            ],
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'vcf_files', to_key => 'vcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3, from_key => 'vcf_files', to_key => 'vcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'filtered_vcf_files', to_key => 'vcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 5, from_key => 'filtered_vcf_files', to_key => 'vcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 6, from_key => 'filtered_vcf_files', to_key => 'vcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'recal_files', to_key => 'recal_files'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'tranches_files', to_key => 'tranches_files'), VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'recalibrated_vcfs', to_key => 'vcf_files'),],
            [VRPipe::StepBehaviourDefiner->new(after_step => 7, behaviour => 'delete_outputs', act_on_steps => [1, 2, 3, 4, 5, 6], regulated_by => 'cleanup', default_regulation => 1),]);
    }
}

1;
