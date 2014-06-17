
=head1 NAME

VRPipe::Steps::convex_gam_correction - a step

=head1 DESCRIPTION

Runs GAMCorrectionPerSample R script from the CoNVex package, generating GAM
Correction files from Read Depth and L2R files

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012-2014 Genome Research Limited.

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

class VRPipe::Steps::convex_gam_correction extends VRPipe::Steps::r_script {
    around options_definition {
        return {
            %{ $self->$orig },
            'convex_rscript_path' => VRPipe::StepOption->create(description => 'full path to CoNVex R scripts'),
            'gc_only'             => VRPipe::StepOption->create(description => 'Boolean for gc_only option in CoNVex gam correction step', default_value => 0),
        };
    }
    
    method inputs_definition {
        return {
            rd_files         => VRPipe::StepIODefinition->create(type => 'rd',  max_files => -1, description => 'read depth (rd) txt files per-sample', metadata => { sample => 'sample name' }),
            l2r_files        => VRPipe::StepIODefinition->create(type => 'l2r', max_files => -1, description => 'l2r txt files per-sample',             metadata => { sample => 'sample name' }),
            features_file    => VRPipe::StepIODefinition->create(type => 'fts', max_files => 1,  description => 'features file from L2R calculation step'),
            breakpoints_file => VRPipe::StepIODefinition->create(type => 'bp',  max_files => 1,  description => 'breakpoints file from convex_breakpoints step'),
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $convex_rscript_path = $options->{'convex_rscript_path'};
            my $gc_only             = $options->{'gc_only'};
            
            my $features_file    = $self->inputs->{features_file}->[0]->path;
            my $breakpoints_file = $self->inputs->{breakpoints_file}->[0]->path;
            
            my %samples;
            foreach my $txt (@{ $self->inputs->{rd_files} }, @{ $self->inputs->{l2r_files} }) {
                my $sample = $txt->metadata->{sample};
                $samples{$sample}{rd}  = $txt       if ($txt->type eq 'rd');
                $samples{$sample}{l2r} = $txt->path if ($txt->type eq 'l2r');
            }
            
            my $req = $self->new_requirements(memory => 2000, time => 1);
            foreach my $sample (keys %samples) {
                my $rd_file  = $samples{$sample}{rd};
                my $l2r_path = $samples{$sample}{l2r};
                
                my $rd_path  = $rd_file->path;
                my $basename = $rd_file->basename;
                $basename =~ s/rd$/gam/;
                
                my $gam_file = $self->output_file(output_key => 'gam_files', basename => $basename, type => 'gam', metadata => $rd_file->metadata);
                my $gam_path = $gam_file->path;
                
                my $cmd = $self->rscript_cmd_prefix . " $convex_rscript_path/GAMCorrectionPerSample.R $l2r_path,$features_file,$gam_path,$rd_path,$breakpoints_file,$gc_only";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::convex_gam_correction', 'run_gam_correction', [$cmd, $req, { output_files => [$gam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return { gam_files => VRPipe::StepIODefinition->create(type => 'gam', max_files => -1, description => 'a GAM Correction file for each input Read depth and L2R file', metadata => { sample => 'sample name' }) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs CoNVex GAMCorrectionPerSample, generating GAM Correction files from Read Depth and L2R files";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method run_gam_correction (ClassName|Object $self: Str $cmd_line) {
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my ($output_path) = $cmd_line =~ / \S+,\S+,(\S+),\S+,\S+$/;
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
        unless ($output_lines > 0) {
            $output_file->unlink;
            $self->throw("Output GAM Correction file is zero length");
        }
        else {
            return 1;
        }
    }
}

1;
