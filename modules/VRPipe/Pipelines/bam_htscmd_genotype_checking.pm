
=head1 NAME

VRPipe::Pipelines::bam_htscmd_genotype_checking - a pipeline

=head1 DESCRIPTION

Uses htscmd gtcheck to check that the genotype of bam files matches the genotype of the samples they claim to be of, 
putting the results of the check into the bam metadata.

=head1 AUTHOR

Chris Joyce <cj5@sangser.ac.uk>.

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

class VRPipe::Pipelines::bam_htscmd_genotype_checking with VRPipe::PipelineRole {
    method name {
        return 'bam_htscmd_genotype_checking';
    }
    
    method description {
        return 'Uses htscmd gtcheck to check that the genotype of bam files matches the genotype of the samples they claim to be of.';
    }
    
    method step_names {
        (
            'bam_index',                 #1
            'vcf_sites',                 #2
            'mpileup_vcf',               #3
            'vcf_index',                 #4
            'htscmd_gtcheck',            #5
            'htscmd_genotype_analysis',  #6
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 0, to_step => 3, to_key   => 'bam_files' },
            { from_step => 1, to_step => 3, from_key => 'bai_files', to_key => 'bai_files' },
            { from_step => 2, to_step => 3, from_key => 'sites_file', to_key   => 'sites_file' },
            { from_step => 3, to_step => 4, from_key => 'vcf_files', to_key => 'vcf_files' },
            { from_step => 3, to_step => 5, from_key => 'vcf_files', to_key => 'vcf_files' },
            { from_step => 5, to_step => 6, from_key => 'htscmd_gtcheck_files', to_key => 'htscmd_gtcheck_files' },
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 6, behaviour => 'delete_outputs', act_on_steps => [2,3,4,5], regulated_by => 'cleanup', default_regulation => 0 });
    }
}

1;
