
=head1 NAME

VRPipe::Pipelines::chipseq_qc_and_peak_calling - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

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

class VRPipe::Pipelines::chipseq_qc_and_peak_calling with VRPipe::PipelineRole {
    method name {
        return 'chipseq_qc_and_peak_calling';
    }
    
    method description {
        return 'Quality check ChIP-seq bams using phantompeakqualtools and macs2.';
    }
    
    method step_names {
        (
            'bam_processing',       #1
            'bam_metadata',         #2
            'phantompeakqualtools', #3
            'macs_callpeak',        #4
            'bedgraph2bigwig',      #5
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 1, to_step => 2, from_key => 'processed_bam_files', to_key => 'bam_files' },
            { from_step => 1, to_step => 3, from_key => 'processed_bam_files', to_key => 'bam_files' },
            { from_step => 1, to_step => 4, from_key => 'processed_bam_files', to_key => 'bam_files' },
            { from_step => 4, to_step => 5, from_key => 'bdg_files', to_key => 'bdg_files' },
        );
    }

}

1;
