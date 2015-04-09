
=head1 NAME

VRPipe::Pipelines::bam_2015_vrtrack_processing_no_recal - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Pipelines::bam_2015_vrtrack_processing_no_recal with VRPipe::PipelineRole {
    method name {
        return 'bam_2015_vrtrack_processing_no_recal';
    }
    
    method description {
        return 'Copy bam files stored in iRODs (found with the irods all_with_warehouse_metadata datasource) to local disc, generate QC stats and graphs, realign around known indels, add NM and BQ tags, checks genotype, does auto QC and updates the VRTrack db. Finally, it archives the bam files. Uses 2015 versions of samtoools, bcftools, picard and GATK.';
    }
    
    method step_names {
        (
            'vrtrack_populate_from_vrpipe_metadata', #1
            'irods_get_files_by_basename',           #2
            'samtools_fasta_gc_stats',               #3
            'samtools_bam_stats',                    #4
            'plot_bamstats',                         #5
            'vrtrack_update_mapstats',               #6
            'sequence_dictionary',                   #7
            'bam_metadata',                          #8 # still necessary?
            'bam_index',                             #9
            'gatk_target_interval_creator',          #10
            'bam_realignment_around_known_indels',   #11
            'bam_calculate_bq',                      #12
            'bam_reheader',                          #13 # still necessary?
            'vrtrack_update_improved',               #14
            'bam_index',                             #15
            'bcftools_generate_sites_file',          #16
            'mpileup_vcf',                           #17
            'vcf_index',                             #18
            'bcftools_gtcheck',                      #19
            'bcftools_genotype_analysis',            #20
            'vrtrack_auto_qc',                       #21
            'archive_files'                          #22
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0,  to_step => 1,  to_key   => 'files' },
            { from_step => 0,  to_step => 2,  to_key   => 'basenames' },
            { from_step => 2,  to_step => 4,  from_key => 'local_files', to_key => 'bam_files' },
            { from_step => 3,  to_step => 5,  from_key => 'fasta_gc_stats_file', to_key => 'fasta_gc_stats_file' },
            { from_step => 4,  to_step => 5,  from_key => 'stats_files', to_key => 'stats_files' },
            { from_step => 2,  to_step => 6,  from_key => 'local_files', to_key => 'bam_files' },
            { from_step => 5,  to_step => 6,  from_key => 'bamstats_plots', to_key => 'bamcheck_plots' },
            { from_step => 2,  to_step => 8,  from_key => 'local_files', to_key => 'bam_files' },
            { from_step => 2,  to_step => 9,  from_key => 'local_files', to_key => 'bam_files' },
            { from_step => 10, to_step => 11, from_key => 'intervals_file', to_key => 'intervals_file' },
            { from_step => 2,  to_step => 11, from_key => 'local_files', to_key => 'bam_files' },
            { from_step => 9,  to_step => 11, from_key => 'bai_files', to_key => 'bai_files' },
            { from_step => 11, to_step => 12, from_key => 'realigned_bam_files', to_key => 'bam_files' },
            { from_step => 12, to_step => 13, from_key => 'bq_bam_files', to_key => 'bam_files' },
            { from_step => 7,  to_step => 13, from_key => 'reference_dict', to_key => 'dict_file' },
            { from_step => 13, to_step => 14, from_key => 'headed_bam_files', to_key => 'bam_files' },
            { from_step => 13, to_step => 15, from_key => 'headed_bam_files', to_key => 'bam_files' },
            { from_step => 13, to_step => 16, from_key => 'headed_bam_files', to_key => 'genotypes_bcf' },
            { from_step => 13, to_step => 17, from_key => 'headed_bam_files', to_key => 'bam_files' },
            { from_step => 15, to_step => 17, from_key => 'bai_files', to_key => 'bai_files' },
            { from_step => 16, to_step => 17, from_key => 'sites_file', to_key => 'sites_file' },
            { from_step => 17, to_step => 18, from_key => 'vcf_files', to_key => 'vcf_files' },
            { from_step => 13, to_step => 19, from_key => 'headed_bam_files', to_key => 'genotypes_bcf' },
            { from_step => 17, to_step => 19, from_key => 'vcf_files', to_key => 'vcf_files' },
            { from_step => 19, to_step => 20, from_key => 'bcftools_gtcheck_files', to_key => 'gtcheck_files' },
            { from_step => 13, to_step => 21, from_key => 'headed_bam_files', to_key => 'bam_files' },
            { from_step => 4,  to_step => 21, from_key => 'stats_files', to_key => 'bamcheck_files' },
            { from_step => 13, to_step => 22, from_key => 'headed_bam_files', to_key => 'file' }
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 11, behaviour => 'delete_outputs', act_on_steps => [2, 9], regulated_by => 'delete_imported_bam', default_regulation => 1 },
            { after_step => 12, behaviour => 'delete_outputs', act_on_steps => [11], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 13, behaviour => 'delete_outputs', act_on_steps => [12], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 21, behaviour => 'delete_outputs', act_on_steps => [17, 18], regulated_by => 'cleanup', default_regulation => 1 }
        );
    }
}

1;
