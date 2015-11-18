
=head1 NAME

VRPipe::Pipelines::vrtrack_auto_qc_from_graph_db - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Pipelines::vrtrack_auto_qc_from_graph_db with VRPipe::PipelineRole {
    method name {
        return 'vrtrack_auto_qc_from_graph_db';
    }
    
    method description {
        return 'Considering the summary stats for a lane in the graph db and the metadata stored on the bam file automatically decide if the lane passes the quality check and update the VRTrack database.';
    }
    
    method step_names {
        (
            'graph_auto_qc',                #1
            'update_vrtrack_from_graph_db', #2
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key => 'bam_files' },
            { from_step => 0, to_step => 2, to_key => 'bam_files' },
        );
    }

}

1;
