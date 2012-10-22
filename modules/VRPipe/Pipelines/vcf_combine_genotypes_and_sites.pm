
=head1 NAME

VRPipe::Pipelines::vcf_combine_genotypes_and_sites - a pipeline

=head1 DESCRIPTION

Merge a sites only VCF file with the genotypes found in another VCF file(s), 
keeping the FILTER, ID and INFO fields found in the sites only file.

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

class VRPipe::Pipelines::vcf_combine_genotypes_and_sites with VRPipe::PipelineRole {
    method name {
        return 'vcf_combine_genotypes_and_sites';
    }
    
    method description {
        return 'Merge a sites only VCF file with the genotypes found in another VCF file(s), keeping the FILTER, ID and INFO fields found in the sites only file';
    }
    
    method step_names {
        (
            'vcf_index',                       #1
            'vcf_combine_genotypes_and_sites', #2
            'vcf_index',                       #3
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'vcf_files' },
            { from_step => 0, to_step => 2, to_key   => 'vcf_files' },
            { from_step => 2, to_step => 3, from_key => 'combined_vcf', to_key => 'vcf_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 3, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'remove_input_vcfs', default_regulation => 0 },
            { after_step => 3, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'remove_input_vcfs', default_regulation => 0 }
        );
    }
}

1;
