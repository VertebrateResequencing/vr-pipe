
=head1 NAME

VRPipe::Pipelines::rna_seq_kallisto_quantification - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2016 Genome Research Limited.

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

class VRPipe::Pipelines::rna_seq_kallisto_quantification with VRPipe::PipelineRole {
    method name {
        return 'rna_seq_kallisto_quantification';
    }
    
    method description {
        return 'Given RNA-seq bam files, quantify transcript abundances with kallisto.';
    }
    
    method step_names {
        (
            'bamtofastq',     #1
            'kallisto_index', #2
            'kallisto_quant', #3
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 1, to_step => 3, from_key => 'fastq_files', to_key => 'fastq_files' },
            { from_step => 2, to_step => 3, from_key => 'transcriptome_index', to_key => 'transcriptome_index' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 3, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1 },
        );
    }
}

1;
