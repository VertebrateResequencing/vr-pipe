
=head1 NAME

VRPipe::Steps::convex_cnv_call - a step

=head1 DESCRIPTION

Runs the SWCNVCall R script in the CoNVex package, generating CNV Calls from
GAM Correction files

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

class VRPipe::Steps::convex_cnv_call extends VRPipe::Steps::r_script {
    around options_definition {
        return {
            %{ $self->$orig },
            'centromere_reg_file' => VRPipe::StepOption->create(description => 'Centromere regions file'),
            'sw_exec'             => VRPipe::StepOption->create(description => 'Full path to Smith-Waterman execution binary'),
            'sw_pval'             => VRPipe::StepOption->create(description => 'p value for Smith-Waterman algorithm', optional => 1, default_value => 2),
            'swt_del'             => VRPipe::StepOption->create(description => 't value threshold for the selection of deletion calls', optional => 1, default_value => 5),
            'swt_dup'             => VRPipe::StepOption->create(description => 't value threshold for the selection of duplication calls', optional => 1, default_value => 5),
            'dv'                  => VRPipe::StepOption->create(description => 'number of probes exponent in CNV call selection', optional => 1, default_value => 0.5),
            'convex_rscript_path' => VRPipe::StepOption->create(description => 'full path to CoNVex R scripts'),
        };
    }
    
    method inputs_definition {
        return { gam_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => '1 or more gam files from which to call CNV') };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $centromere_reg_file = $options->{'centromere_reg_file'};
            my $sw_exec             = $options->{'sw_exec'};
            my $sw_pval             = $options->{'sw_pval'};
            my $swt_del             = $options->{'swt_del'};
            my $swt_dup             = $options->{'swt_dup'};
            my $dv                  = $options->{'dv'};
            my $convex_rscript_path = $options->{'convex_rscript_path'};
            
            my $req = $self->new_requirements(memory => 2000, time => 1);
            
            foreach my $gam_file (@{ $self->inputs->{gam_files} }) {
                my $gam_path = $gam_file->path;
                my $sample   = $gam_file->metadata->{sample};
                
                my $basename = $gam_file->basename;
                $basename =~ s/gam\.txt$/cnv_calls.txt/;
                
                my $cnv_file = $self->output_file(output_key => 'cnv_files', basename => $basename, type => 'txt');
                my $cnv_dir = $cnv_file->dir;
                
                my $cmd = $self->rscript_cmd_prefix . " $convex_rscript_path/SWCNVCall.R $sw_pval,$swt_del,$swt_dup,$dv,$gam_path,$basename,$sample,$centromere_reg_file,$cnv_dir,$sw_exec";
                $self->dispatch([$cmd, $req, { output_files => [$cnv_file] }]);
            }
        };
    
    }
    
    method outputs_definition {
        return { cnv_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'a CNV call file for each input GAM file'), };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs CoNVex SW CNV Call program, generating a CNV Call text file for each GAM file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }

}

1;
