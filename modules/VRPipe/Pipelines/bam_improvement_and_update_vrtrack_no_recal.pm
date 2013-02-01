
=head1 NAME

VRPipe::Pipelines::bam_improvement_and_update_vrtrack_no_recal - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012-2013 Genome Research Limited.

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

class VRPipe::Pipelines::bam_improvement_and_update_vrtrack_no_recal with VRPipe::PipelineRole {
    method name {
        return 'bam_improvement_and_update_vrtrack_no_recal';
    }
    
    method description {
        return 'Improves bam files by realigning around known indels, and adding NM and BQ tags; adds the resulting bam to your VRTrack db. Does no quality score recalibration, useful when dealing with species that diverge from reference greatly or have no/few SNP calls';
    }
    
    method step_names {
        (
            'sequence_dictionary',                 #1
            'bam_metadata',                        #2
            'bam_index',                           #3
            'gatk_target_interval_creator',        #4
            'bam_realignment_around_known_indels', #5
            'bam_calculate_bq',                    #6
            'bam_reheader',                        #7
            'vrtrack_update_improved',             #8
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 2, to_key   => 'bam_files' },
            { from_step => 0, to_step => 3, to_key   => 'bam_files' },
            { from_step => 4, to_step => 5, from_key => 'intervals_file', to_key => 'intervals_file' },
            { from_step => 0, to_step => 5, to_key   => 'bam_files' },
            { from_step => 3, to_step => 5, from_key => 'bai_files', to_key => 'bai_files' },
            { from_step => 5, to_step => 6, from_key => 'realigned_bam_files', to_key => 'bam_files' },
            { from_step => 6, to_step => 7, from_key => 'bq_bam_files', to_key => 'bam_files' },
            { from_step => 1, to_step => 7, from_key => 'reference_dict', to_key => 'dict_file' },
            { from_step => 7, to_step => 8, from_key => 'headed_bam_files', to_key => 'bam_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 5, behaviour => 'delete_inputs',  act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0 },
            { after_step => 5, behaviour => 'delete_outputs', act_on_steps => [3], regulated_by => 'delete_input_bams', default_regulation => 0 },
            { after_step => 6, behaviour => 'delete_outputs', act_on_steps => [4, 5], regulated_by => 'cleanup', default_regulation => 1 },
            { after_step => 7, behaviour => 'delete_outputs', act_on_steps => [6], regulated_by => 'cleanup', default_regulation => 1 }
        );
    }
}

1;
