
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

class VRPipe::Steps::convex_plots extends VRPipe::Steps::r_script {
    method inputs_definition {
        return { cnv_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'Set of convex CNV call files '), };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my @o_files;
            push @o_files, $self->output_file(output_key => 'calls_per_sample_png', basename => "CNVstats_CallsperSample.png", type => 'bin');
            push @o_files, $self->output_file(output_key => 'del_dup_ratio_png',    basename => "CNVstats_DelDupRatio.png",    type => 'bin');
            
            # Create a fofn of CNV call files
            my $cnv_fofn = $self->output_file(basename => "cnv.fofn", type => 'txt', temporary => 1);
            $cnv_fofn->create_fofn($self->inputs->{cnv_files});
            
            my $cnv_fofn_path = $cnv_fofn->path;
            my $out_dir       = $cnv_fofn->dir;
            
            my $plot_script_path = $self->output_file(basename => "convex_plot.R", type => 'txt', temporary => 1)->path;
            
            my $req = $self->new_requirements(memory => 1200, time => 1);
            my $cmd = $self->rscript_cmd_prefix . " $plot_script_path $cnv_fofn_path,$out_dir";
            
            $self->dispatch_wrapped_cmd('VRPipe::Steps::convex_plots', 'run_plots', [$cmd, $req, { output_files => \@o_files }]);
        };
    }
    
    method outputs_definition {
        return {
            calls_per_sample_png => VRPipe::StepIODefinition->create(type => 'bin', max_files => 1, description => 'CNVstats_CallsperSample png'),
            del_dup_ratio_png    => VRPipe::StepIODefinition->create(type => 'bin', max_files => 1, description => 'CNVstats_DelDupRatio png'),
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
    
    method run_plots (ClassName|Object $self: Str $cmd_line!) {
        my ($rscript_path) = $cmd_line =~ m/(\S+\.R) \S+$/;
        $rscript_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $rscript = VRPipe::File->get(path => $rscript_path);
        
        my $fh = $rscript->openw;
        print $fh qq[
require(CoNVex)
ca = commandArgs(trailingOnly=TRUE);

fn = strsplit(ca,"\\\\,")
fofn = fn[[1]][1]
outdir = fn[[1]][2]

setwd(outdir)

CNVfiles = as.character(read.table(fofn)[,1])
CNVcalls = GetCNVCalls(CNVfiles)

PlotCNVStats(CNVcallsAll=CNVcalls)
];
        $rscript->close;
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        return 1;
    }

}

1;
