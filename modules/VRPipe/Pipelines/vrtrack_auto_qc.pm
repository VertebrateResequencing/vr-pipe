=head1 NAME

VRPipe::Pipelines::vrtrack_auto_qc - a pipeline

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

class VRPipe::Pipelines::vrtrack_auto_qc with VRPipe::PipelineRole {
    method name {
        return 'vrtrack_auto_qc';
    }
    method _num_steps {
        return 1;
    }
    method description {
        return 'Considering the stats in the bamcheck file for a lane, and the metadata stored on the bam file and in the VRTrack database for the corresponding lane, automatically decide if the lane passes the quality check.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'vrtrack_auto_qc') ],
                [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                  VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bamcheck_files') ],
                [ ]);
    }
}

1;