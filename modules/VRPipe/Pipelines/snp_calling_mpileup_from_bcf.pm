
=head1 NAME

VRPipe::Pipelines::§ - a pipeline

=head1 DESCRIPTION

Runs samtools mpileup for a bam datasource, generating both bcf and (using
bcftools view) vcf files. Run this if you need to keep the intermediate bcfs,
otherwise use the snp_calling_mpileup_vcf pipeline.

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Pipelines::snp_calling_mpileup_from_bcf with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_mpileup_from_bcf';
    }
    
    method description {
        return 'Run bcftools to produce vcf from bcf files';
    }
    
    method step_names {
        (
            'bcf_to_vcf',
            'vcf_index'
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bcf_files' },
            { from_step => 0, to_step => 1, to_key   => 'sites_file' },
            { from_step => 1, to_step => 2, from_key => 'vcf_files', to_key => 'vcf_files' }
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 2, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'remove_bcfs', default_regulation => 0 });
    }
}

1;
