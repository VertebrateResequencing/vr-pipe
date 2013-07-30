
=head1 NAME

VRPipe::Pipelines::cnv_control_comparison - a pipeline

=head1 DESCRIPTION

Remove control CNV from stem cell CNVs.

=head1 AUTHOR

Phil Carter <pc12@sanger.ac.uk>, John Maslen <jm23@sanger.ac.uk>.

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

class VRPipe::Pipelines::cnv_control_comparison with VRPipe::PipelineRole {
    method name {
        return 'cnv_control_comparison';
    }
    
    method description {
        return 'Reformat CNV files from penncnv or quantisnp into bed format, then compare sample CNV output to that of the control cell type. Produces output files of intersection between control and sample and sample minus contol (i.e. the diff)';
    }
    
    method step_names {
        (
            'reformat_cnv_output_to_bed',
            'cnv_control_comparison'
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'cnv_files' },
            { from_step => 1, to_step => 2, from_key => 'bed_files', to_key => 'bed_files' }
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 0 });
    }
}

1;
