
=head1 NAME

VRPipe::Pipelines::bam_merge_and_streaming_mark_duplicates - a pipeline

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

class VRPipe::Pipelines::bam_merge_and_streaming_mark_duplicates with VRPipe::PipelineRole {
    method name {
        return 'bam_merge_and_streaming_mark_duplicates';
    }
    
    method description {
        return 'Merge and mark duplicates with biobambam bamstreamingmarkduplicates. Input files must have the tags required by bamstreamingmarkduplicates';
    }
    
    method step_names {
        (
            'samtools_merge_and_streaming_mark_duplicates', #1
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key => 'bam_files' },
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 1, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0 });
    }
}

1;
