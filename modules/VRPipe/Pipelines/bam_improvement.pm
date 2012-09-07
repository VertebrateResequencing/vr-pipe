
=head1 NAME

VRPipe::Pipelines::bam_improvement - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Pipelines::bam_improvement with VRPipe::PipelineRole {
    method name {
        return 'bam_improvement';
    }
    
    method _num_steps {
        return 10;
    }
    
    method description {
        return 'Improves bam files by realigning around known indels, recalibrating quality scores and adding NM and BQ tags';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
                VRPipe::Step->get(name => 'sequence_dictionary'),                 #1
                VRPipe::Step->get(name => 'bam_metadata'),                        #2
                VRPipe::Step->get(name => 'bam_index'),                           #3
                VRPipe::Step->get(name => 'gatk_target_interval_creator'),        #4
                VRPipe::Step->get(name => 'bam_realignment_around_known_indels'), #5
                VRPipe::Step->get(name => 'bam_index'),                           #6
                VRPipe::Step->get(name => 'bam_count_covariates'),                #7
                VRPipe::Step->get(name => 'bam_recalibrate_quality_scores'),      #8
                VRPipe::Step->get(name => 'bam_calculate_bq'),                    #9
                VRPipe::Step->get(name => 'bam_reheader'),                        #10
            ],
            
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'intervals_file', to_key => 'intervals_file'), VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 5, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'realigned_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 7, from_key => 'realigned_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 8, from_key => 'bam_recalibration_files', to_key => 'bam_recalibration_files'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 8, from_key => 'realigned_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 8, to_step => 9, from_key => 'recalibrated_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 9, to_step => 10, from_key => 'bq_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 10, from_key => 'reference_dict', to_key => 'dict_file')],
            
            [VRPipe::StepBehaviourDefiner->new(after_step => 5, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0), VRPipe::StepBehaviourDefiner->new(after_step => 5, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'delete_input_bams', default_regulation => 0), VRPipe::StepBehaviourDefiner->new(after_step => 8, behaviour => 'delete_outputs', act_on_steps => [5, 6, 7], regulated_by => 'cleanup', default_regulation => 1), VRPipe::StepBehaviourDefiner->new(after_step => 9, behaviour => 'delete_outputs', act_on_steps => [8], regulated_by => 'cleanup', default_regulation => 1), VRPipe::StepBehaviourDefiner->new(after_step => 10, behaviour => 'delete_outputs', act_on_steps => [9], regulated_by => 'cleanup', default_regulation => 1)]
        );
    }
}

1;
