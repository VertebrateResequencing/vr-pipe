
=head1 NAME

VRPipe::Pipelines::bam_improvement_no_recal - a pipeline

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

class VRPipe::Pipelines::bam_improvement_no_recal with VRPipe::PipelineRole {
    method name {
        return 'bam_improvement_no_recal';
    }
    
    method _num_steps {
        return 7;
    }
    
    method description {
        return 'Improves bam files by realigning around known indels, and adding NM and BQ tags. Does no quality score recalibration, useful when dealing with species that diverge from reference greatly or have no/few SNP calls';
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
             VRPipe::Step->get(name => 'bam_calculate_bq'),                    #6
             VRPipe::Step->get(name => 'bam_reheader'),                        #7
            ],
            
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'intervals_file', to_key => 'intervals_file'), VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 5, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'realigned_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'bq_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 7, from_key => 'reference_dict', to_key => 'dict_file')],
            
            [VRPipe::StepBehaviourDefiner->new(after_step => 5, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0), VRPipe::StepBehaviourDefiner->new(after_step => 5, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'delete_input_bams', default_regulation => 0), VRPipe::StepBehaviourDefiner->new(after_step => 6, behaviour => 'delete_outputs', act_on_steps => [4, 5], regulated_by => 'cleanup', default_regulation => 1), VRPipe::StepBehaviourDefiner->new(after_step => 7, behaviour => 'delete_outputs', act_on_steps => [6], regulated_by => 'cleanup', default_regulation => 1)]);
    }
}

1;
