
=head1 NAME

VRPipe::Pipelines::bam_lanelet_gt_check - a pipeline

=head1 DESCRIPTION

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

class VRPipe::Pipelines::bam_lanelet_gt_check with VRPipe::PipelineRole {
    method name {
        return 'bam_lanelet_gt_check';
    }
    
    method description {
        return 'Runs lanelet Genotype consistency check and updates bam metadata if GT is confirmed';
    }
    
    method step_names {
        (
            'lanelet_gt_bam_select',     #1
            'bam_index',                 #2
            'vcf_sites',                 #3
            'mpileup_vcf',               #4
            'vcf_index',                 #5
            'htscmd_gtcheck',            #6
            'lanelet_gt_bam_update',     #7
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'bam_files' },
            { from_step => 1, to_step => 2, from_key => 'bam_files', to_key  => 'bam_files' },
            { from_step => 1, to_step => 4, from_key => 'bam_files', to_key  => 'bam_files' },
            { from_step => 2, to_step => 4, from_key => 'bai_files', to_key => 'bai_files' },
            { from_step => 3, to_step => 4, from_key => 'sites_file', to_key => 'sites_file' },
            { from_step => 4, to_step => 5, from_key => 'vcf_files', to_key => 'vcf_files' },
            { from_step => 4, to_step => 6, from_key => 'vcf_files', to_key => 'vcf_files' },
            { from_step => 6, to_step => 7, from_key => 'htscmd_gtcheck_files', to_key => 'htscmd_gtcheck_files' },
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 7, behaviour => 'delete_outputs', act_on_steps => [1,2,3,4,5,6], regulated_by => 'cleanup', default_regulation => 1 });
    }
}

1;
