
=head1 NAME

VRPipe::Pipelines::bam_genotype_checking - a pipeline

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

class VRPipe::Pipelines::bam_genotype_checking with VRPipe::PipelineRole {
    method name {
        return 'bam_genotype_checking';
    }
    
    method _num_steps {
        return 4;
    }
    
    method description {
        return 'Check that the genotype of bam files matches the genotype of the samples they claim to be of.';
    }
    
    method steps {
        $self->throw("steps cannot be called on this non-persistent object");
    }
    
    method _step_list {
        return ([
                VRPipe::Step->get(name => 'bin2hapmap_sites'),         #1
                VRPipe::Step->get(name => 'mpileup_bcf_hapmap'),       #2
                VRPipe::Step->get(name => 'glf_check_genotype'),       #3
                VRPipe::Step->get(name => 'gtypex_genotype_analysis'), #4
            ],
            
            [VRPipe::StepAdaptorDefiner->new(from_step => 0, to_step => 2, to_key => 'bam_files'), VRPipe::StepAdaptorDefiner->new(from_step => 1, to_step => 2, from_key => 'hapmap_file', to_key => 'hapmap_file'), VRPipe::StepAdaptorDefiner->new(from_step => 2, to_step => 3, from_key => 'bcf_files_with_metadata', to_key => 'bcf_files'), VRPipe::StepAdaptorDefiner->new(from_step => 3, to_step => 4, from_key => 'gtypex_files_with_metadata', to_key => 'gtypex_files'),],
            
            [VRPipe::StepBehaviourDefiner->new(after_step => 4, behaviour => 'delete_outputs', act_on_steps => [1, 2], regulated_by => 'cleanup', default_regulation => 0),]
        );
    }
}

1;
