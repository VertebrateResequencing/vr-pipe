
=head1 NAME

VRPipe::Pipelines::bam_import_from_irods - a pipeline

=head1 DESCRIPTION

One step pipeline to retrieve bams from irods.

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

class VRPipe::Pipelines::bam_import_from_irods with VRPipe::PipelineRole {
    method name {
        return 'bam_import_from_irods';
    }
    
    method _num_steps {
        return 1;
    }
    
    method description {
        return 'Copy bam files stored in iRODs to local disc.';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return (
            [VRPipe::Step->get(name => 'irods_get_files_by_basename')],
            
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'basenames')],
            
            []
        );
    }
}

1;
