
=head1 NAME

VRPipe::Pipelines::convex_plot_generation - a pipeline

=head1 DESCRIPTION

This pipeline generates CNV call plots from Convex CNV call text files. The 
pipeline runs once only for a set of cnv calls; its datasource will probably be
a vrpipe group_by_metadata set of text files from a convex_cnv_calling pipeline
or your own fofn of cnv calls.

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

class VRPipe::Pipelines::convex_plot_generation with VRPipe::PipelineRole {
    method name {
        return 'convex_plot_generation';
    }
    
    method description {
        return 'Run CoNVex pipeline to Generate CNV call plots from set of Convex CNV calls';
    }
    
    method step_names {
        ('convex_plots');
    }
    
    method adaptor_definitions {
        ({ from_step => 0, to_step => 1, to_key => 'cnv_files' });
    }
}

1;
