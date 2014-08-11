
=head1 NAME

VRPipe::Pipelines::shapeit - a pipeline

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

class VRPipe::Pipelines::shapeit with VRPipe::PipelineRole {
    method name {
        return 'shapeit';
    }
    
    method description {
        return 'Run SHAPEIT to prephase a VCF file for imputation.';
    }
    
    method step_names {
        (
            'define_vcf_chunks',           #1
            'create_imputation_ref_panel', #2
            'shapeit',                     #3
            'bcftools_concat',             #4
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'vcf_files' },
            { from_step => 0, to_step => 3, to_key   => 'vcf_files' },
            { from_step => 1, to_step => 2, from_key => 'bed_files', to_key => 'bed_files' },
            { from_step => 1, to_step => 3, from_key => 'bed_files', to_key => 'bed_files' },
            { from_step => 3, to_step => 4, from_key => 'vcf_files', to_key => 'vcf_files' },
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 4, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'cleanup', default_regulation => 0 });
    }

}

1;
