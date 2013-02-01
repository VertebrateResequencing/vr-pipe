
=head1 NAME

VRPipe::Pipelines::snp_calling_mpileup - a pipeline

=head1 DESCRIPTION

Runs samtools mpileup followed by bcftools view, generating vcf files for a bam
datasource. Run the snp_calling_mpileup_via_bcf pipeline instead if you need to
keep the intermediate bcfs.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Pipelines::snp_calling_mpileup with VRPipe::PipelineRole {
    method name {
        return 'snp_calling_mpileup';
    }
    
    method description {
        return 'Run samtools mpileup to generate vcf files';
    }
    
    method step_names {
        (
            'bam_index',   #1
            'mpileup_vcf', #2
            'vcf_index'    #3
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 0, to_step => 2, to_key   => 'bam_files' },
            { from_step => 1, to_step => 2, from_key => 'bai_files', to_key => 'bai_files' },
            { from_step => 2, to_step => 3, from_key => 'vcf_files', to_key => 'vcf_files' }
        );
    }
}

1;
