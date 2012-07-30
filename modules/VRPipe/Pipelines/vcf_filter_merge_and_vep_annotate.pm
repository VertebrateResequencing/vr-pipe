
=head1 NAME

VRPipe::Pipelines::vcf_filter_merge_and_vep_annotate - a pipeline

=head1 DESCRIPTION

Pipeline to filter, merge and annotate VCF files with consequences. The datasource will be a delimited datsource, each element consting of two bam files. These are initially filtered indepenedently, then merged using vcf-isec, and annotated using vcf-annotate. The vcfs are then annotated with consequences using the Ensembl Variant Effect Predictor, VEP, and stats are generated for the final vcfs.

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

class VRPipe::Pipelines::vcf_filter_merge_and_vep_annotate with VRPipe::PipelineRole {
    method name {
        return 'vcf_filter_merge_and_vep_annotate';
    }
    
    method _num_steps {
        return 6;
    }
    
    method description {
        return 'Filter then merge VCF files and then annotate them with consequences using VEP';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
             VRPipe::Step->get(name => 'vcf_multi_filter'),     #1
             VRPipe::Step->get(name => 'vcf_merge'),            #2
             VRPipe::Step->get(name => 'vcf_annotate'),         #3
             VRPipe::Step->get(name => 'vep_analysis'),         #4
             VRPipe::Step->get(name => 'vcf_vep_consequences'), #5
             VRPipe::Step->get(name => 'vcf_stats')],                                                  #6
            
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'vcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'filtered_vcf', to_key => 'vcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'merged_vcf', to_key => 'vcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'annotated_vcf', to_key => 'vcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 5, from_key => 'annotated_vcf', to_key => 'vcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'vep_txt', to_key => 'vep_txt'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'conseq_vcf', to_key => 'vcf_files')],
            
            [VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1), VRPipe::StepBehaviourDefiner->new(after_step => 3, behaviour => 'delete_outputs', act_on_steps => [2], regulated_by => 'cleanup', default_regulation => 1), VRPipe::StepBehaviourDefiner->new(after_step => 5, behaviour => 'delete_outputs', act_on_steps => [3, 4], regulated_by => 'cleanup', default_regulation => 1)]);
    }
}

1;
