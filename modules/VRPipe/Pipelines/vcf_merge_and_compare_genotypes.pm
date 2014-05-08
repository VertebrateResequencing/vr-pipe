
=head1 NAME

VRPipe::Pipelines::vcf_merge_and_compare_genotypes - a pipeline

=head1 DESCRIPTION

Pipeline to merge multiple VCF files with vcftools vcf-isec, then compare the
genotypes of the samples to see if they all came from the same individual.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Pipelines::vcf_merge_and_compare_genotypes with VRPipe::PipelineRole {
    method name {
        return 'vcf_merge_and_compare_genotypes';
    }
    
    method description {
        return 'Merge multiple VCF files with vcf-isec, then compare genotypes of the samples to see if they all came from the same individual';
    }
    
    method step_names {
        (
            'vcf_merge_different_samples', #1
            'vcf_genotype_comparison',     #2
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'vcf_files' },
            { from_step => 1, to_step => 2, from_key => 'merged_vcf', to_key => 'vcf_files' },
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 2, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'remove_input_vcfs', default_regulation => 0 });
    }
}

1;
