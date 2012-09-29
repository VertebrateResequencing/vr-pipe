
=head1 NAME

VRPipe::Pipelines::vcf_filter_merge_and_vep_annotate - a pipeline

=head1 DESCRIPTION

Pipeline to filter, merge and annotate VCF files with consequences. The
datasource will be a delimited datsource, each element consting of two bam
files. These are initially filtered indepenedently, then merged using vcf-isec,
and annotated using vcf-annotate. The vcfs are then annotated with consequences
using the Ensembl Variant Effect Predictor, VEP, and stats are generated for
the final vcfs.

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

class VRPipe::Pipelines::vcf_filter_merge_and_vep_annotate with VRPipe::PipelineRole {
    method name {
        return 'vcf_filter_merge_and_vep_annotate';
    }
    
    method description {
        return 'Filter then merge VCF files and then annotate them with consequences using VEP';
    }
    
    method step_names {
        (
            'vcf_multi_filter',     #1
            'vcf_index',            #2
            'vcf_merge',            #3
            'vcf_index',            #4
            'vcf_annotate',         #5
            'vcf_index',            #6
            'vep_analysis',         #7
            'vcf_vep_consequences', #8
            'vcf_index',            #9
            'vcf_stats'             #10
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1,  to_key   => 'vcf_files' },
            { from_step => 1, to_step => 2,  from_key => 'filtered_vcf', to_key => 'vcf_files' },
            { from_step => 1, to_step => 3,  from_key => 'filtered_vcf', to_key => 'vcf_files' },
            { from_step => 3, to_step => 4,  from_key => 'merged_vcf', to_key => 'vcf_files' },
            { from_step => 3, to_step => 5,  from_key => 'merged_vcf', to_key => 'vcf_files' },
            { from_step => 5, to_step => 6,  from_key => 'annotated_vcf', to_key => 'vcf_files' },
            { from_step => 5, to_step => 7,  from_key => 'annotated_vcf', to_key => 'vcf_files' },
            { from_step => 5, to_step => 8,  from_key => 'annotated_vcf', to_key => 'vcf_files' },
            { from_step => 7, to_step => 8,  from_key => 'vep_txt', to_key => 'vep_txt' },
            { from_step => 8, to_step => 9,  from_key => 'conseq_vcf', to_key => 'vcf_files' },
            { from_step => 8, to_step => 10, from_key => 'conseq_vcf', to_key => 'vcf_files' }
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 2, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'remove_input_vcfs', default_regulation => 0 },
            { after_step => 4, behaviour => 'delete_outputs', act_on_steps => [1, 2], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 6, behaviour => 'delete_outputs', act_on_steps => [3, 4], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 9, behaviour => 'delete_outputs', act_on_steps => [5, 6, 7], regulated_by => 'cleanup', default_regulation => 1 }
        );
    }
}

1;
