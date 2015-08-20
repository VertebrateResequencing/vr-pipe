
=head1 NAME

VRPipe::Pipelines::bcf_processing - a pipeline

=head1 DESCRIPTION

One step pipeline for the ad-hoc bcf_processing step.

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

class VRPipe::Pipelines::bcf_processing with VRPipe::PipelineRole {
    method name {
        return 'bcf_processing';
    }
    
    method description {
        return 'Run a VCF or BCF file through some ad-hoc command line to produce another VCF or BCF file';
    }
    
    method step_names {
        (
            'bcf_processing', #1
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key => 'bcf_files' },
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 1, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'delete_input_bcfs', default_regulation => 0 });
    }
}

1;
