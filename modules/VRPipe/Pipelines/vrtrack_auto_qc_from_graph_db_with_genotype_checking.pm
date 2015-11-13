
=head1 NAME

VRPipe::Pipelines::vrtrack_auto_qc_from_graph_db_with_genotype_checking - a
pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Pipelines::vrtrack_auto_qc_from_graph_db_with_genotype_checking with VRPipe::PipelineRole {
    method name {
        return 'vrtrack_auto_qc_from_graph_db_with_genotype_checking';
    }
    
    method description {
        return 'Uses bcftools gtcheck to check that the genotype of bam files matches the genotype of the samples they claim to be of. Then, considering the summary stats for a lane in the graph db and the metadata stored on the bam file automatically decide if the lane passes the quality check and update the VRTrack database.';
    }
    
    method step_names {
        (
            
            'bcftools_generate_sites_file', #1
            'mpileup_vcf',                  #2
            'bcftools_gtcheck',             #3
            'bcftools_genotype_analysis',   #4
            'graph_auto_qc',                #5
            'update_vrtrack_from_graph_db', #6
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'genotypes_bcf' },
            { from_step => 0, to_step => 2, to_key   => 'bam_files' },
            { from_step => 0, to_step => 2, to_key   => 'bai_files' },
            { from_step => 1, to_step => 2, from_key => 'sites_file', to_key => 'sites_file' },
            { from_step => 0, to_step => 3, to_key   => 'genotypes_bcf' },
            { from_step => 2, to_step => 3, from_key => 'vcf_files', to_key => 'vcf_files' },
            { from_step => 3, to_step => 4, from_key => 'bcftools_gtcheck_files', to_key => 'gtcheck_files' },
            { from_step => 0, to_step => 5, to_key   => 'bam_files' },
            { from_step => 0, to_step => 6, to_key   => 'bam_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 4, behaviour => 'delete_outputs', act_on_steps => [2, 3], regulated_by => 'cleanup', default_regulation => 1 },
        );
    }

}

1;
