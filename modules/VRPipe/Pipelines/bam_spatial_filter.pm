
=head1 NAME

VRPipe::Pipelines::bam_spatial_filter - a pipeline

=head1 DESCRIPTION



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

class VRPipe::Pipelines::bam_spatial_filter with VRPipe::PipelineRole {
    method name {
        return 'bam_spatial_filter';
    }
    
    method description {
        return 'apply spatial filter to bams';
    }
    
    method step_names {
        (
            'calculate_bam_spatial_filter',
            'apply_bam_spatial_filter',
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step  => 1,              to_key  => 'bam_files' },
            { from_step => 1, from_key => 'filter_files', to_step => 2, to_key => 'filter_files' },
            { from_step => 0, to_step  => 2,              to_key  => 'bam_files' },
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1 });
    }
}

1;
