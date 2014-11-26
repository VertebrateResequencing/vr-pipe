
=head1 NAME

VRPipe::Pipelines::pluritest_gene_expression_analysis - a pipeline

=head1 DESCRIPTION

Reformats the gene expression files outputted by Genome Studio and then runs
the R PluriTest  analysis package on these reformatted files to generate graphs
and ancillary data for analysis of pluripotency.

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>. Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013,2014 Genome Research Limited.

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

class VRPipe::Pipelines::pluritest_gene_expression_analysis with VRPipe::PipelineRole {
    method name {
        return 'pluritest_gene_expression_analysis';
    }
    
    method description {
        return 'Reformat the Genome Studio csv files for gene expression data and run the R PluriTest analysis package to generate graphs and ancillary data';
    }
    
    method step_names {
        (
            'pluritest_annotation_profile_files',
            'pluritest_reformat_genome_studio_expression_files',
            'pluritest_plot_gene_expression'
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'idat_files' },
            { from_step => 1, to_step => 2, from_key => 'profile_file', to_key => 'profile_files' },
            { from_step => 1, to_step => 2, from_key => 'annotation_file', to_key => 'annotation_files' },
            { from_step => 1, to_step => 2, from_key => 'mapping_file', to_key => 'mapping_files' },
            { from_step => 2, to_step => 3, from_key => 'reformat_files', to_key => 'conv_files' }
        );
    }
    
    method behaviour_definitions {
        ({ after_step => 3, behaviour => 'delete_outputs', act_on_steps => [1, 2], regulated_by => 'cleanup', default_regulation => 0 });
    }
}

1;
