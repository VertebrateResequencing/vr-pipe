
=head1 NAME

VRPipe::Pipelines::convex_mean_mad_calculation - a pipeline

=head1 DESCRIPTION



=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Pipelines::convex_mean_mad_calculation with VRPipe::PipelineRole {
    method name {
        return 'convex_mean_mad_calculation';
    }
    
    method description {
        return 'Add mean and mad statistics to convex CNV calls and generate convex plots';
    }
    
    method step_names {
        (
            'convex_mean_mad', #1
            'convex_plots',    #2
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 0, to_step => 1, to_key   => 'rd_files' },
            { from_step => 0, to_step => 1, to_key   => 'l2r_files' },
            { from_step => 0, to_step => 1, to_key   => 'gam_files' },
            { from_step => 0, to_step => 1, to_key   => 'cnv_files' },
            { from_step => 1, to_step => 2, from_key => 'cnv_files_with_mean_mad', to_key => 'cnv_files' }
        );
    }
}

1;
