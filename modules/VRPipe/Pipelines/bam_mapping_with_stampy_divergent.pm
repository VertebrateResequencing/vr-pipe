=head1 NAME

VRPipe::Pipelines::bam_mapping_with_stampy_divergent - a pipeline

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

class VRPipe::Pipelines::bam_mapping_with_stampy_divergent with VRPipe::PipelineRole {
    method name {
        return 'bam_mapping_with_stampy_divergent';
    }
    method _num_steps {
        return 15;
    }
    method description {
        return 'Map reads in bam files to a reference genome with stampy, using a preliminary mapping and snp calling to determine the substitution rate to map with';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'sequence_dictionary'), #1
                  VRPipe::Step->get(name => 'stampy_buildgenome'), #2
                  VRPipe::Step->get(name => 'stampy_buildhash'), #3
                  VRPipe::Step->get(name => 'bam_metadata'), #4
                  VRPipe::Step->get(name => 'bam_name_sort'), #5
                  VRPipe::Step->get(name => 'bam_to_fastq'), #6
                  VRPipe::Step->get(name => 'fastq_split'), #7
                  VRPipe::Step->get(name => 'stampy_map_fastq'), #8
                  VRPipe::Step->get(name => 'sam_to_fixed_bam'), #9
                  VRPipe::Step->get(name => 'bam_merge_lane_splits'), #10
                  VRPipe::Step->get(name => 'bam_substitution_rate'), #11
                  VRPipe::Step->get(name => 'stampy_map_fastq'), #12
                  VRPipe::Step->get(name => 'sam_to_fixed_bam'), #13
                  VRPipe::Step->get(name => 'bam_merge_lane_splits'), #14
                  VRPipe::Step->get(name => 'bamcheck') #15
                 ],
                 
                 [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 4, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 5, to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'name_sorted_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 8, from_key => 'split_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 8, to_step => 9, from_key => 'stampy_sam_files', to_key => 'sam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 9, to_step => 10, from_key => 'fixed_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 10, from_key => 'reference_dict', to_key => 'dict_file'),
                   
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 11, from_key => 'reference_dict', to_key => 'dict_file'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 10, to_step => 11, from_key => 'merged_lane_bams', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 12, from_key => 'split_fastq_files', to_key => 'fastq_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 12, to_step => 13, from_key => 'stampy_sam_files', to_key => 'sam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 13, to_step => 14, from_key => 'fixed_bam_files', to_key => 'bam_files'),
                   VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 14, from_key => 'reference_dict', to_key => 'dict_file'),
                   
                   VRPipe::StepAdaptorDefiner->new(from_step => 14, to_step => 15, from_key => 'merged_lane_bams', to_key => 'bam_files') ],
                 
                 [ VRPipe::StepBehaviourDefiner->new(after_step => 6, behaviour => 'delete_outputs', act_on_steps => [5], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 9, behaviour => 'delete_outputs', act_on_steps => [8], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 10, behaviour => 'delete_outputs', act_on_steps => [9], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 11, behaviour => 'delete_outputs', act_on_steps => [10], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 12, behaviour => 'delete_outputs', act_on_steps => [6, 7], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 13, behaviour => 'delete_outputs', act_on_steps => [12], regulated_by => 'cleanup', default_regulation => 1),
                   VRPipe::StepBehaviourDefiner->new(after_step => 14, behaviour => 'delete_outputs', act_on_steps => [13], regulated_by => 'cleanup', default_regulation => 1) ]);
    }
}

1;
