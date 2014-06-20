
=head1 NAME

VRPipe::Steps::verify_bamid - a step

=head1 DESCRIPTION

Runs verifyBamID as a contamination check to verify whether the reads in BAM
files match the genotypes in SNP-only VCF files.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>, Shane McCarthy <sm15@sanger.ac.uk>

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

class VRPipe::Steps::verify_bamid with VRPipe::StepRole {
    method options_definition {
        return {
            verify_bamid_exe  => VRPipe::StepOption->create(description => 'path to verifyBamID executable', optional => 1, default_value => 'verifyBamID'),
            verify_bamid_opts => VRPipe::StepOption->create(description => 'verifyBamID options excluding --bam and --out. Must inlcude the --vcf option'),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                description => '1 or more bam files',
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
            $self->throw("Invalid options: do not include --bam or --out options. '$verify_bamid_opts'") if ($verify_bamid_opts =~ /--bam/ or $verify_bamid_opts =~ /--out/);
            $self->throw("Invalid options: must include the --vcf option. '$verify_bamid_opts'") unless ($verify_bamid_opts =~ /--vcf/);
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'verifyBamID',
                    version => VRPipe::StepCmdSummary->determine_version($verify_bamid_exe, '^verifyBamID (\S+)'),
                    summary => 'verifyBamID ' . $verify_bamid_opts . ' --bam $bam --out $out_prefix'
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $bam_file (@{ $self->inputs->{bam_files} }) {
                my $bam_path = $bam_file->path;
                my $basename = $bam_file->basename;
                $basename =~ s/\.bam$//;
                
                my $out_log = $self->output_file(output_key => 'verify_bam_id_logs', basename => "$basename.log", type => 'txt', metadata => { source_bam_id => $bam_file->id });
                my $out_path = $out_log->path;
                $out_path =~ s/\.log$//;
                
                my $self_sm_file  = $self->output_file(output_key => 'verify_bam_id_self_sm_files',  basename => "$basename.selfSM",  type => 'txt', metadata => { source_bam_id => $bam_file->id });
                my $depth_sm_file = $self->output_file(output_key => 'verify_bam_id_depth_sm_files', basename => "$basename.depthSM", type => 'txt', metadata => { source_bam_id => $bam_file->id });
                my @output_files = ($out_log, $self_sm_file, $depth_sm_file);
                
                unless ($verify_bamid_opts =~ /ignoreRG/) {
                    my $self_rg_file = $self->output_file(output_key => 'verify_bam_id_self_rg_files', basename => "$basename.selfRG", type => 'txt', metadata => { source_bam_id => $bam_file->id });
                    push(@output_files, $self_rg_file);
                    my $depth_rg_file = $self->output_file(output_key => 'verify_bam_id_depth_rg_files', basename => "$basename.depthRG", type => 'txt', metadata => { source_bam_id => $bam_file->id });
                    push(@output_files, $depth_rg_file);
                }
                if ($verify_bamid_opts =~ /best/) {
                    my $best_sm_file = $self->output_file(output_key => 'verify_bam_id_best_sm_files', basename => "$basename.bestSM", type => 'txt', metadata => { source_bam_id => $bam_file->id });
                    push(@output_files, $best_sm_file);
                    unless ($verify_bamid_opts =~ /ignoreRG/) {
                        my $best_rg_file = $self->output_file(output_key => 'verify_bam_id_best_rg_files', basename => "$basename.bestRG", type => 'txt', metadata => { source_bam_id => $bam_file->id });
                        push(@output_files, $best_rg_file);
                    }
                }
                
                my $cmd = qq[$verify_bamid_exe $verify_bamid_opts --bam $bam_path --out $out_path];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::verify_bamid', 'verify_bam_id_and_check', [$cmd, $req, { output_files => \@output_files }]);
            }
        
        };
    }
    
    method outputs_definition {
        return {
            verify_bam_id_self_sm_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID sample match stats',
                max_files   => -1,
                min_files   => 0
            ),
            verify_bam_id_depth_sm_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID depth distribution stats',
                max_files   => -1
            ),
            verify_bam_id_logs => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID log file',
                max_files   => -1
            ),
            verify_bam_id_self_rg_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID per-readgroup lane-sample match stats',
                max_files   => -1,
                min_files   => 0
            ),
            verify_bam_id_depth_rg_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID per-readgroup depth distribution stats',
                max_files   => -1,
                min_files   => 0
            ),
            verify_bam_id_best_sm_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID sample best-match stats',
                max_files   => -1,
                min_files   => 0
            ),
            verify_bam_id_best_rg_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'verifyBamID per-readgroup sample best-match stats',
                max_files   => -1,
                min_files   => 0
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs verifyBamID contamination check on input BAM files; Adds results like FREEMIX and CHIPMIX from selfSM stats file as metadata to the BAM files.";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method verify_bam_id_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($bam_path, $out_path) = $cmd_line =~ /--bam (\S+) --out (\S+)$/;
        $bam_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $log_file = VRPipe::File->get(path => "$out_path.log");
        my $selfSM   = VRPipe::File->get(path => "$out_path.selfSM");
        my $bam_file = VRPipe::File->get(path => $bam_path);
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $log_file->update_stats_from_disc;
        my $logh = $log_file->openr;
        my $ok   = 0;
        while (my $line = <$logh>) {
            if ($line =~ /^Analysis finished/) {
                $ok++;
            }
        }
        $logh->close;
        $self->throw("Analysis not completed, see $out_path.log") unless $ok;
        
        # #SEQ_ID RG  CHIP_ID #SNPS   #READS  AVG_DP  FREEMIX FREELK1 FREELK0 FREE_RH FREE_RA CHIPMIX CHIPLK1 CHIPLK0 CHIP_RH CHIP_RA DPREF   RDPHET  RDPALT
        # SC_AutoBirm5533219  ALL NA  81610   2286086 28.01   0.00298 511343.41   513182.69   0.51918 0.00003 NA  NA  NA  NA  NA  NA  NA  NA
        $selfSM->update_stats_from_disc;
        my $fh = $selfSM->openr;
        my @cols;
        my %meta;
        while (<$fh>) {
            chomp;
            $self->throw("Unexpected number of rows in verifyBamID output file. selfSM files expected to have only one sample: " . $selfSM->path) if (keys %meta);
            if (/^#SEQ_ID/) {
                @cols = split /\t/;
            }
            else {
                $self->throw("Unexpected number of columns in verifyBamID output file: " . $selfSM->path) unless (@cols == 19);
                my @vals = split /\t/;
                foreach my $col (@cols) {
                    my $val = shift @vals;
                    next if ($col eq '#SEQ_ID' || $col eq 'RG');
                    next if $val eq 'NA';
                    $col =~ s/#/N_/;
                    $meta{"VBID_$col"} = $val;
                }
            }
        }
        $fh->close;
        
        $bam_file->add_metadata(\%meta);
        
        return 1;
    }

}

1;
