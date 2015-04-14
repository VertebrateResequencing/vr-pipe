
=head1 NAME

VRPipe::Pipelines::vcf_concat - a pipeline

=head1 DESCRIPTION

Pipeline to concatenate multiple vcf files using the bcftools concat.

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

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

class VRPipe::Pipelines::bcftools_concat with VRPipe::PipelineRole {
    method name {
        return 'bcftools_concat';
    }
    
    method description {
        return 'Concatenate multiple vcfs using bcftools concat';
    }
    
    method step_names {
        (
            'bcftools_concat', #1
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key => 'vcf_files' },
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 1, behaviour => 'delete_inputs', act_on_steps => [0], regulated_by => 'remove_input_vcfs', default_regulation => 0 });
    }
}

1;
