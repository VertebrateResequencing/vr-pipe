
=head1 NAME

VRPipe::Pipelines::genome_studio_import_from_irods_and_convert_to_vcf - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <jm23@sanger.ac.uk>.

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

class VRPipe::Pipelines::genome_studio_import_from_irods_and_convert_to_vcf with VRPipe::PipelineRole {
    method name {
        return 'genome_studio_import_from_irods_and_convert_to_vcf';
    }
    
    method description {
        return 'Import gtc files from iRODS along with metadata containing location of stored genome studio genotypes file. The metadata is used to obtain individual sample genotype files for genotype analysis. Also converts the genotype files (fcr) to VCF.';
    }
    
    method step_names {
        (
            'irods_get_files_by_basename',        #1
            'split_genome_studio_genotype_files', #2
            'illumina_coreexome_manifest_to_map', #3
            'genome_studio_fcr_to_vcf',           #4
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'basenames' },
            { from_step => 1, to_step => 2, from_key => 'local_files', to_key => 'gtc_files' },
            { from_step => 2, to_step => 4, from_key => 'gtype_files', to_key => 'fcr_files' },
            { from_step => 3, to_step => 4, from_key => 'map_file', to_key => 'map_file' }
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1 },
        );
    }

}

1;
