
=head1 NAME

VRPipe::Pipelines::fastq_mapping_with_smalt - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Pipelines::fastq_mapping_with_smalt with VRPipe::PipelineRole {
    method name {
        return 'fastq_mapping_with_smalt';
    }
    
    method _num_steps {
        return 10;
    }
    
    method description {
        return 'Map reads in fastq files to a reference genome with smalt';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
             VRPipe::Step->get(name => 'sequence_dictionary'),   #1
             VRPipe::Step->get(name => 'smalt_index'),           #2
             VRPipe::Step->get(name => 'fastq_import'),          #3
             VRPipe::Step->get(name => 'fastq_metadata'),        #4
             VRPipe::Step->get(name => 'fastq_split'),           #5
             VRPipe::Step->get(name => 'fastq_decompress'),      #6
             VRPipe::Step->get(name => 'smalt_map_to_sam'),      #7
             VRPipe::Step->get(name => 'sam_to_fixed_bam'),      #8
             VRPipe::Step->get(name => 'bam_add_readgroup'),     #9
             VRPipe::Step->get(name => 'bam_merge_lane_splits'), #10
            ],
            
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'fastq_files'), VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'local_fastq_files', to_key => 'fastq_files'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'fastq_files_with_metadata', to_key => 'fastq_files'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'split_fastq_files', to_key => 'compressed_fastq_files'), VRPipe::StepAdaptorDefiner->new(from_step => 6, to_step => 7, from_key => 'decompressed_fastq_files', to_key => 'fastq_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 7, from_key => 'reference_dict', to_key => 'dict_file'), VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 7, from_key => 'smalt_index_binary_files', to_key => 'index_files'), VRPipe::StepAdaptorDefiner->new(from_step => 7, to_step => 8, from_key => 'smalt_sam_files', to_key => 'sam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 8, to_step => 9, from_key => 'fixed_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 9, to_step => 10, from_key => 'rg_added_bam_files', to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 10, from_key => 'reference_dict', to_key => 'dict_file')],
            
            [VRPipe::StepBehaviourDefiner->new(after_step => 7, behaviour => 'delete_outputs', act_on_steps => [3, 5, 6], regulated_by => 'cleanup', default_regulation => 1), VRPipe::StepBehaviourDefiner->new(after_step => 8, behaviour => 'delete_outputs', act_on_steps => [7], regulated_by => 'cleanup', default_regulation => 1), VRPipe::StepBehaviourDefiner->new(after_step => 9, behaviour => 'delete_outputs', act_on_steps => [8], regulated_by => 'cleanup', default_regulation => 1), VRPipe::StepBehaviourDefiner->new(after_step => 10, behaviour => 'delete_outputs', act_on_steps => [9], regulated_by => 'cleanup', default_regulation => 1)]);
    }
}

1;
