
=head1 NAME

VRPipe::Pipelines::sga_variant_calling - a pipeline

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

class VRPipe::Pipelines::sga_variant_calling with VRPipe::PipelineRole {
    method name {
        return 'sga_variant_calling';
    }
    
    method _num_steps {
        return 6;
    }
    
    method description {
        return 'Use sga to call variants. Merge and index fastq, sga index and sga graph-diff.';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
                VRPipe::Step->get(name => 'sga_permute_reference'),             #1
                VRPipe::Step->get(name => 'sga_index_reference'),               #2
                VRPipe::Step->get(name => 'sga_generate_sampled_suffix_array'), #3
                VRPipe::Step->get(name => 'fastq_merge_and_index'),             #4
                VRPipe::Step->get(name => 'sga_index'),                         #5
                VRPipe::Step->get(name => 'sga_reference_based_calling'),       #6
            ],
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 4, to_key => 'fastq_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'permuted_reference_fasta', to_key => 'reference_fasta'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3, from_key => 'permuted_reference_fasta', to_key => 'reference_fasta'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'merged_fastq_file', to_key => 'fastq_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 6, from_key => 'permuted_reference_fasta', to_key => 'reference_fasta'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 6, from_key => 'merged_fastq_file', to_key => 'sga_indexed_variant_reads'),],
            [VRPipe::StepBehaviourDefiner->new(after_step => 4, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'remove_inputs', default_regulation => 0), VRPipe::StepBehaviourDefiner->new(after_step => 6, behaviour => 'delete_outputs', act_on_steps => [4, 5], regulated_by => 'cleanup', default_regulation => 0),]
        );
    }
}

1;
