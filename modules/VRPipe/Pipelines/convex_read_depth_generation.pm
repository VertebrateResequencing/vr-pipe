
=head1 NAME

VRPipe::Pipelines::convex_read_depth_generation - a pipeline

=head1 DESCRIPTION

This is first of three pipeline required to run in sequence in order to
generate CNV Calls using the Convex Exome CNV detection package. This pipeline
generates Read Depth files from a BAM Datasource, and adds Sample and Sex
metadata to the bams, for subsequent L2R Generation and CNV Calling pipelines.

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

class VRPipe::Pipelines::convex_read_depth_generation with VRPipe::PipelineRole {
    method name {
        return 'convex_read_depth_generation';
    }
    
    method description {
        return 'Run CoNVex pipeline to Generate Read Depth files and add Sample metadata to bam files, for subsequent L2R Generation and CNV Calling pipelines';
    }
    
    method step_names {
        (
            'convex_read_depth',     #1
            'bam_metadata_with_sex', #2
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key => 'bam_files' },
            { from_step => 0, to_step => 2, to_key => 'bam_files' }
        );
    }
}

1;
