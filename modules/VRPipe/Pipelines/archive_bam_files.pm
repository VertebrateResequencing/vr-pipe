
=head1 NAME

VRPipe::Pipelines::archive_bam_files - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Pipelines::archive_bam_files with VRPipe::PipelineRole {
    method name {
        return 'archive_bam_files';
    }
    
    method description {
        return 'Safely move BAM/CRAM and BAI/CRAI index files from one disc to a pool of one or more other discs, eg. for archival purposes. (the output_root option for this pipeline is meaningless)';
    }
    
    method step_names {
        (
            'archive_bam_files',       #1
            'archive_bam_index_files', #2
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key => 'bam_file' },
            { from_step => 0, to_step => 2, to_key => 'bam_file' },
            { from_step => 0, to_step => 2, to_key => 'bam_index_file' },
        );
    }
}

1;
