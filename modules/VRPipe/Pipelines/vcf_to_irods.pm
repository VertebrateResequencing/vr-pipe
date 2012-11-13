
=head1 NAME

VRPipe::Pipelines::vcf_to_irods - a pipeline

=head1 DESCRIPTION

Pipeline to add vr-pipe metadata to vcf files, and then add them into irods with the required metadata.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Pipelines::vcf_to_irods with VRPipe::PipelineRole {
    method name {
        return 'vcf_to_irods';
    }
    
    method description {
        return 'Add vcf files with metadata to irods';
    }
    
    method step_names {
        (
            'vcf_metadata', #1
            'vcf_to_irods', #2
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key => 'vcf_files' },
            { from_step => 0, to_step => 2, to_key => 'vcf_files' },

        );
    }
    
}

1;
