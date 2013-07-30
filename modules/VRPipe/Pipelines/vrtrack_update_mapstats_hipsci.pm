
=head1 NAME

VRPipe::Pipelines::vrtrack_update_mapstats_hipsci - a pipeline

=head1 DESCRIPTION

Update vrtrack mapstats table with cnv calling data. Specifically, numbers of
CNVs with and without those  intersecting with the control sample CNVs removed.
The chromosome coordinates are also stored.

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Pipelines::vrtrack_update_mapstats_hipsci with VRPipe::PipelineRole {
    method name {
        return 'vrtrack_update_mapstats_hipsci';
    }
    
    method description {
        return 'Pipeline to update vrtrack mapstats for hipsci genotyping post CNV calling';
    }
    
    method step_names {
        ('vrtrack_update_mapstats_hipsci');
    }
    
    method adaptor_definitions {
        ({ from_step => 0, to_step => 1, to_key => 'diff_files' });
    }
}

1;
