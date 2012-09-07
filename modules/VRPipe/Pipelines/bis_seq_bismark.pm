
=head1 NAME

VRPipe::Pipelines::bis-seq_bismark - a pipeline

=head1 DESCRIPTION

Bisulfite Sequencing pipeline using Bismark and associated software
(http://www.bioinformatics.babraham.ac.uk/projects/bismark/) to produce
methylation calls.

=head1 AUTHOR

NJWalker <nw11@sanger.ac.uk>.

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

class VRPipe::Pipelines::bis_seq_bismark with VRPipe::PipelineRole {
    method name {
        return 'bis_seq_bismark';
    }
    
    method _num_steps {
        return 4;
    }
    
    method description {
        return 'Bisulphite Sequencing Pipeline employing Bismark tools to gain methylation calls.';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
                VRPipe::Step->get(name => 'fastqc_quality_report'),        # 1
                VRPipe::Step->get(name => 'trimmomatic'),                  # 2
                VRPipe::Step->get(name => 'bismark'),                      # 3
                VRPipe::Step->get(name => 'sam_sort'),                     # 4
                VRPipe::Step->get(name => 'sam_mark_duplicates'),          # 5
                VRPipe::Step->get(name => 'bismark_methylation_extractor') # 6
            ],
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'fastq_files'), VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'fastq_files'), VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'trimmed_files', to_key => 'fastq_files'), VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'bismark_sam', to_key => 'sam_file'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'sorted_sam', to_key => 'sam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 4, to_step => 5, from_key => 'sorted_sam', to_key => 'sam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 5, to_step => 6, from_key => 'markdup_sam_files', to_key => 'sam_file')],
            [
            
            ]
        );
    }
}
1;
