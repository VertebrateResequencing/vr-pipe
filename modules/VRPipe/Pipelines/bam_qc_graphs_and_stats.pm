=head1 NAME

VRPipe::Pipelines::bam_qc_graphs_and_stats - a pipeline

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

class VRPipe::Pipelines::bam_qc_graphs_and_stats with VRPipe::PipelineRole {
    method name {
        return 'bam_qc_graphs_and_stats';
    }
    method _num_steps {
        return 3;
    }
    method description {
        return 'Generates graphs and stats to describe bam files for quality-checking purposes.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'fasta_gc_stats'),#1
                  VRPipe::Step->get(name => 'bamcheck'),#2
                  VRPipe::Step->get(name => 'plot_bamcheck'),#3
                ],
                
                [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'),
                  VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 3, from_key => 'fasta_gc_stats_file', to_key => 'fasta_gc_stats_file'),
                  VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'bamcheck_files', to_key => 'bamcheck_files')
                ],
                [ ]
               );
    }
}

1;