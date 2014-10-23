
=head1 NAME

VRPipe::Pipelines::irods_analysis_files_download - a pipeline

=head1 DESCRIPTION

Pipeline to download all files from an irods all_with_warehouse_metadata
datasource.

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

class VRPipe::Pipelines::irods_analysis_files_download with VRPipe::PipelineRole {
    method name {
        return 'irods_analysis_files_download';
    }
    
    method description {
        return 'Download files from irods specified by an irods all_with_warehouse_metadata datasource - gets associated analysis files and (optionally) the files listed by your query.';
    }
    
    method step_names {
        ('irods_analysis_files_download');
    }
    
    method adaptor_definitions {
        ({ from_step => 0, to_step => 1, to_key => 'files' });
    }
}

1;
