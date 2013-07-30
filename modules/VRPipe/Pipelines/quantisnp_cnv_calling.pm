
=head1 NAME

VRPipe::Pipelines::quantisnp_cnv_calling - a pipeline

=head1 DESCRIPTION

Reformats the gene expression files outputted by Genome Studio and then runs
the R PluriTest  analysis package on these reformatted files to generate graphs
and ancillary data for analysis of pluripotency.

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

class VRPipe::Pipelines::quantisnp_cnv_calling with VRPipe::PipelineRole {
    method name {
        return 'quantisnp_cnv_calling';
    }
    
    method description {
        return 'Reformat genotyping data for use with QuantiSNP, then detect CNVs using Quantisnp';
    }
    
    method step_names {
        (
            'quantisnp_reformat_gs_export',
            'quantisnp_detect_cnv'
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key => 'stepOne_file_input_GS_file' },                                                              # 1st step takes the gs file as input
            { from_step => 1, to_step => 2, from_key =>, 'stepOne_file_output_reformatted_file', to_key => 'stepTwo_file_input_reformatted_file' } # 2nd step takes as input the reformat files produced by the first step
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 0 });
    }
}

1;                                                                                                                                                 # needs this to exit correctly as a package
