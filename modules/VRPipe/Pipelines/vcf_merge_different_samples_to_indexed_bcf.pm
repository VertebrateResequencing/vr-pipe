
=head1 NAME

VRPipe::Pipelines::vcf_merge - a pipeline

=head1 DESCRIPTION

Pipeline to merge multiple VCF files with different samples to an indexed BCF

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Pipelines::vcf_merge_different_samples_to_indexed_bcf with VRPipe::PipelineRole {
    method name {
        return 'vcf_merge_different_samples_to_indexed_bcf';
    }
    
    method description {
        return 'Merge multiple VCF files for different samples to an indexed BCF file';
    }
    
    method step_names {
        (
            'vcf_merge_different_samples_to_indexed_bcf' #1
        );
    }
    
    method adaptor_definitions {
        ({ from_step => 0, to_step => 1, to_key => 'vcf_files' });
    }
}

1;
