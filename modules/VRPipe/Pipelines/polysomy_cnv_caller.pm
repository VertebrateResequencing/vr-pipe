
=head1 NAME

VRPipe::Pipelines::polysomy_cnv_caller - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Pipelines::polysomy_cnv_caller with VRPipe::PipelineRole {
    method name {
        return 'polysomy_cnv_caller';
    }
    
    method description {
        return 'Merges different sample VCFs together and detects trisomy/tetrasomy and copy number variants; requires Illumina B-allele frequency scores (BAF) and a control sample in input VCF files.';
    }
    
    method step_names {
        (
            'vcf_merge_different_samples_control_aware', #1
            'polysomy',                                  #2
            'plot_polysomy',                             #3
            'bcftools_cnv',                              #4
            'combine_bcftools_cnvs'                      #5
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'vcf_files' },
            { from_step => 1, to_step => 2, from_key => 'merged_vcf', to_key => 'vcf_files' },
            { from_step => 1, to_step => 4, from_key => 'merged_vcf', to_key => 'vcf_files' },
            { from_step => 2, to_step => 3, from_key => 'dist_file', to_key => 'dist_files' },
            { from_step => 4, to_step => 5, from_key => 'summary_file', to_key => 'summary_files' },
        );
    }
}

1;
