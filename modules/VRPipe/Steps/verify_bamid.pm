
=head1 NAME

VRPipe::Steps::verify_bamid - a step

=head1 DESCRIPTION

Runs verifyBamID as a contamination check to verify whether the reads in bam
files match the genotypes in SNP-only vcf files Requires that the vcf file path
be a metadata item for each bam

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

class VRPipe::Steps::verify_bamid with VRPipe::StepRole {
    method options_definition {
        return {
            verify_bamid_exe  => VRPipe::StepOption->create(description => 'path to verifyBamID executable',                      optional => 1, default_value => 'verifyBamID'),
            verify_bamid_opts => VRPipe::StepOption->create(description => 'verifyBamID options excluding --vcf --bam and --out', optional => 1, default_value => '--ignoreRG'),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                description => '1 or more bam files',
                metadata    => { vcf => 'path to vcf' },
                max_files   => -1,
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self              = shift;
            my $options           = $self->options;
            my $verify_bamid_exe  = $options->{verify_bamid_exe};
            my $verify_bamid_opts = $options->{verify_bamid_opts};
            $self->throw("Invalid options '$verify_bamid_opts'") if $verify_bamid_opts =~ /--vcf/ or $verify_bamid_opts =~ /--bam/ or $verify_bamid_opts =~ /--out/;
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $bam_file (@{ $self->inputs->{bam_files} }) {
                my $vcf_path = $bam_file->metadata->{vcf};
                my $bam_path = $bam_file->path;
                my $basename = $bam_file->basename;
                $basename =~ s/\.bam$/.vfb/;
                
                my $out_log = $self->output_file(output_key => 'out_logs', basename => "$basename.log", type => 'txt');
                my $out_path = $out_log->path;
                $out_path =~ s/\.log$//;
                
                my $self_sm_file  = $self->output_file(output_key => 'self_sm_files',  basename => "$basename.selfSM",  type => 'txt');
                my $depth_sm_file = $self->output_file(output_key => 'depth_sm_files', basename => "$basename.depthSM", type => 'txt');
                my @output_files = ($out_log, $self_sm_file, $depth_sm_file);
                
                unless ($verify_bamid_opts =~ /ignoreRG/) {
                    my $self_rg_file = $self->output_file(output_key => 'self_rg_files', basename => "$basename.selfRG", type => 'txt');
                    push(@output_files, $self_rg_file);
                    my $depth_rg_file = $self->output_file(output_key => 'depth_rg_files', basename => "$basename.depthRG", type => 'txt');
                    push(@output_files, $depth_rg_file);
                }
                if ($verify_bamid_opts =~ /best/) {
                    my $best_sm_file = $self->output_file(output_key => 'best_sm_files', basename => "$basename.bestSM", type => 'txt');
                    push(@output_files, $best_sm_file);
                    unless ($verify_bamid_opts =~ /ignoreRG/) {
                        my $best_rg_file = $self->output_file(output_key => 'best_rg_files', basename => "$basename.bestRG", type => 'txt');
                        push(@output_files, $best_rg_file);
                    }
                }
                
                my $cmd = qq[use VRPipe::Steps::verify_bamid; VRPipe::Steps::verify_bamid->verify_bam(verify_bamid_exe => '$verify_bamid_exe', verify_bamid_opts => '$verify_bamid_opts', bam_path => '$bam_path', vcf_path => '$vcf_path', out_path => '$out_path');];
                
                $self->dispatch_vrpipecode($cmd, $req, { output_files => \@output_files });
            }
        
        };
    }
    
    method outputs_definition {
        return {
            self_sm_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID sample match stats',
                max_files   => -1,
                min_files   => 0
            ),
            depth_sm_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID depth distribution stats',
                max_files   => -1
            ),
            out_logs => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID  log',
                max_files   => -1
            ),
            self_rg_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID per-readgroup lane-sample match stats',
                max_files   => -1,
                min_files   => 0
            ),
            depth_rg_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID per-readgroup depth distribution stats',
                max_files   => -1,
                min_files   => 0
            ),
            best_sm_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID sample best-match stats',
                max_files   => -1,
                min_files   => 0
            ),
            best_rg_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID er-readgroup sample best-match stats',
                max_files   => -1,
                min_files   => 0
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs verifyBamID contamination check to verify genotype in bam files";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method verify_bam (ClassName|Object $self: Str :$verify_bamid_exe!, Str :$verify_bamid_opts!, Str|File :$bam_path!,  Str|File :$vcf_path!,  Str|File :$out_path! ) {
        my $cmd_line = "$verify_bamid_exe --vcf $vcf_path --bam $bam_path --out $out_path $verify_bamid_opts";
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        $self->warn($cmd_line);
        
        my $log_file = VRPipe::File->get(path => "$out_path.log");
        
        $log_file->update_stats_from_disc;
        my $logh = $log_file->openr;
        my $ok   = 0;
        while (my $line = <$logh>) {
            if ($line =~ /^Analysis finished/) {
                $ok++;
            }
        }
        $self->throw("Analysis not completed, see $out_path.log") unless $ok;
        return 1;
    }

}

1;
