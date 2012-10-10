
=head1 NAME

VRPipe::Steps::convex_plots - a step

=head1 DESCRIPTION

Runs the PlotCNVStats R script from the CoNVex packages. Runs once per
pipeline, generating a set of png plots for the called cnvs.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Steps::convex_plots with VRPipe::StepRole {
    method options_definition {
        return {
            'r_libs' => VRPipe::StepOption->create(description => 'full path to R_LIBS where the CoNVex package is installed'),
        };
    }
    
    method inputs_definition {
        return { cnv_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'Set of convex CNV call files '), };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            my $r_libs  = $options->{'r_libs'};
            
            $self->output_file(output_key => 'calls_per_sample_png', basename => "CNVstats_CallsperSample.png", type => 'bin');
            $self->output_file(output_key => 'del_dup_ratio_png',    basename => "CNVstats_DelDupRatio.png",    type => 'bin');
            
            # Create a fofn of CNV call files
            my $cnv_fofn = $self->output_file(output_key => 'cnv_fofn', basename => "CnvFofn.txt", type => 'txt');
            my $ofh = $cnv_fofn->openw;
            
            foreach my $cnv_file (@{ $self->inputs->{cnv_files} }) {
                my $rd_path = $cnv_file->path;
                print $ofh "$rd_path\n";
            }
            $cnv_fofn->close;
            
            eval "use Statistics::R;";
            my $R = Statistics::R->new();
            $R->run(".libPaths('$r_libs')");
            $R->run("require(CoNVex)");
            
            my $cnv_fofn_path = $cnv_fofn->path;
            $R->run("fofn='$cnv_fofn_path'");
            $R->run("CNVfiles = as.character(read.table(fofn)[,1])");
            $R->run("CNVcalls = GetCNVCalls(CNVfiles)");
            
            my $out_dir = $cnv_fofn->dir;
            $R->run("setwd('$out_dir')");
            $R->run("PlotCNVStats(CNVcallsAll=CNVcalls)");
        };
    }
    
    method outputs_definition {
        return {
            cnv_fofn             => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'fofn of convex call files'),
            calls_per_sample_png => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'CNVstats_CallsperSample png'),
            del_dup_ratio_png    => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'CNVstats_DelDupRatio png'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs CoNVex PlotCNVStats, generating png plots from a set of cnv calls files";
    }
    
    method max_simultaneous {
        return 1;            # should only run once
    }
}

1;
