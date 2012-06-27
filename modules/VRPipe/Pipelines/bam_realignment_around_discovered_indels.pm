=head1 NAME

VRPipe::Pipelines::bam_realignment_around_discovered_indels - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Pipelines::bam_realignment_around_discovered_indels with VRPipe::PipelineRole {
    method name {
        return 'bam_realignment_around_discovered_indels';
    }
    method _num_steps {
        return 4;
    }
    method description {
        return 'Use GATK to discover possible indels in a bam, then realign in those regions.';
    }
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([ VRPipe::Step->get(name => 'bam_metadata'),#1
                  VRPipe::Step->get(name => 'bam_index'),#2
                  VRPipe::Step->get(name => 'gatk_target_interval_creator_discovery'),#3
		  VRPipe::Step->get(name => 'bam_realignment_around_discovered_indels'),#4
		],
		
   	        [ VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 1, to_key => 'bam_files'),
                  VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'),
		  VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 3, to_key => 'bam_files'),
		  VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 4, to_key => 'bam_files'),
		  VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'intervals_file', to_key => 'intervals_file'),
		],
		
		[ VRPipe::StepBehaviourDefiner->new(after_step => 2, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_input_bams', default_regulation => 0),
		  VRPipe::StepBehaviourDefiner->new(after_step => 4, behaviour => 'delete_outputs', act_on_steps => [2, 3], regulated_by => 'cleanup', default_regulation => 1)]);
    }
}

1;
