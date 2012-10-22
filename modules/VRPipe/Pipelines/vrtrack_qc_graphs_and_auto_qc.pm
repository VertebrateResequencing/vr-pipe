
=head1 NAME

VRPipe::Pipelines::vrtrack_qc_graphs_and_auto_qc - a pipeline

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

class VRPipe::Pipelines::vrtrack_qc_graphs_and_auto_qc with VRPipe::PipelineRole {
    method name {
        return 'vrtrack_qc_graphs_and_auto_qc';
    }
    
    method description {
        return 'Given bam files already on disc, generate QC stats and graphs, do the auto-qc anaylsis on the results, and update the VRTrack db.';
    }
    
    method step_names {
        (
            'fasta_gc_stats',          #1
            'bamcheck',                #2
            'plot_bamcheck',           #3
            'vrtrack_update_mapstats', #4
            'vrtrack_auto_qc'          #5
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 2, to_key   => 'bam_files' },
            { from_step => 1, to_step => 3, from_key => 'fasta_gc_stats_file', to_key => 'fasta_gc_stats_file' },
            { from_step => 2, to_step => 3, from_key => 'bamcheck_files', to_key => 'bamcheck_files' },
            { from_step => 0, to_step => 4, to_key   => 'bam_files' },
            { from_step => 3, to_step => 4, from_key => 'bamcheck_plots', to_key => 'bamcheck_plots' },
            { from_step => 0, to_step => 5, to_key   => 'bam_files' },
            { from_step => 2, to_step => 5, from_key => 'bamcheck_files', to_key => 'bamcheck_files' },
        );
    }
}

1;
