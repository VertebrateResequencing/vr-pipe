
=head1 NAME

VRPipe::Pipelines::vrtrack_populate_from_irods_and_download_files - a pipeline

=head1 DESCRIPTION

Pipeline to get metadata from an irods all_with_warehouse_metadata datasource
and populate/update a VRTrack database. It also downloads associated files from
irods, and optionally downloads the primary files as well.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Pipelines::vrtrack_populate_from_irods_and_download_files with VRPipe::PipelineRole {
    method name {
        return 'vrtrack_populate_from_irods_and_download_files';
    }
    
    method description {
        return 'Populate or update VRTrack lanes based on input file metadata from the irods all_with_warehouse_metadata datasource, then download the files from irods.';
    }
    
    method step_names {
        ('vrtrack_populate_from_vrpipe_metadata', 'irods_analysis_files_download');
    }
    
    method adaptor_definitions {
        ({ from_step => 0, to_step => 1, to_key => 'files' }, { from_step => 0, to_step => 2, to_key => 'files' });
    }
}

1;
