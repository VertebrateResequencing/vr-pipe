
=head1 NAME

VRPipe::Steps::macs_callpeak - a step

=head1 DESCRIPTION

This step implements macs2 callpeak for ChIP-seq peak-finding. Input files are 
pre-filtered chipseq bams grouped by individual. E.g. input bams may include 
ChIP data for sample coxy_33 produced using four antibodies:

path     individual      sample file1.bam    coxy_33   coxy_33_K27AC file2.bam 
  coxy_33   coxy_33_K27ME3 file3.bam    coxy_33   coxy_33_K4ME3 file4.bam   
coxy_33   coxy_33_INPUT

Peaks are then called using: macs2 callpeak -t file1.bam -c file4.bam $options1
-n coxy_33_K27AC macs2 callpeak -t file2.bam -c file4.bam $options2 -n
coxy_33_K27ME3 macs2 callpeak -t file3.bam -c file4.bam $options3 -n
coxy_33_K4ME3 macs2 callpeak -t file4.bam -c file4.bam $options4 -n
coxy_33_INPUT

Where there is no control (_INPUT) peak-calling is run without -c.

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::macs_callpeak extends VRPipe::Steps::r {
    around options_definition {
        return {
            %{ $self->$orig },
            macs2_exe              => VRPipe::StepOption->create(description => 'path to your MACS2 executable',                                                                                                                         optional => 1, default_value => 'macs2'),
            sample_metadata_key    => VRPipe::StepOption->create(description => 'sample metadata key that determines the antibody used for enrichment (e.g. ChIPseq of H3K4me3 markers may include the string H3K4me3 in its metadata)', optional => 1),
            control_sample_regex   => VRPipe::StepOption->create(description => 'regular expression string that determines the control sample',                                                                                          optional => 1, default_value => 'INPUT'),
            macs2_callpeak_options => VRPipe::StepOption->create(description => 'options to macs2 callpeak (a comma separated list of marker:macs2_options would be acceptable when peak calling depends on the antibody used)',         optional => 1, default_value => 'K27AC:-g hs,K4ME3:-g hs,K27ME3:-g hs --broad,INPUT:-g hs --broad'),
            bedtools_exe           => VRPipe::StepOption->create(description => 'path to your bedtools executable',                                                                                                                      optional => 1, default_value => 'bedtools'),
            samtools_exe           => VRPipe::StepOption->create(description => 'path to your samtools executable',                                                                                                                      optional => 1, default_value => 'samtools'),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'pre-filtered bams',
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            my $macs2_exe     = $options->{macs2_exe};
            my $macs2_opts    = $options->{macs2_callpeak_options};
            my $sample_key    = $options->{sample_metadata_key};
            my $control_regex = $options->{control_sample_regex};
            my $bedtools_exe  = $options->{bedtools_exe};
            my $samtools_exe  = $options->{samtools_exe};
            
            if ($macs2_opts =~ /\s-[tcn]\s/) {
                $self->throw("macs2_callpeak_options should not include -t, -c or -n");
            }
            my %marker_options;
            my @macs2_opts = split(/,/, $macs2_opts);
            foreach my $opts (@macs2_opts) {
                my @opts = split(/:/, $opts);
                $marker_options{ $opts[0] } = $opts[1];
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => "$macs2_exe",
                    version => VRPipe::StepCmdSummary->determine_version($macs2_exe, '^macs2 (.+)$'),
                    summary => "macs2 callpeak -B -t \$ChIP.bam -c \$Control.bam -n \$sample_name \$macs2_options"
                )
            );
            
            my @input_bams = @{ $self->inputs->{bam_files} };
            my %input_bams;
            foreach my $bam_file (@input_bams) {
                my $sample = $bam_file->metadata->{$sample_key};
                if (exists $input_bams{$sample}) {
                    $self->throw("two bam files with the same sample $sample metadata: $input_bams{$sample}\t$bam_file not allowed!");
                }
                else {
                    $input_bams{$sample} = $bam_file;
                }
            }
            
            my $control_bam;
            foreach my $sample (keys %input_bams) {
                if ($sample =~ /$control_regex/i) {
                    $control_bam = $input_bams{$sample}->path;
                }
            }
            
            my $r_cmd_prefix = $self->r_cmd_prefix;
            foreach my $input_bam (@{ $self->inputs->{bam_files} }) {
                my $sample = $input_bam->metadata->{$sample_key};
                
                my $xls_file        = $self->output_file(output_key => 'xls_file',  basename => "${sample}_peaks.xls",        type => 'txt', metadata => { $sample_key => $sample });
                my $narrowPeak_file = $self->output_file(output_key => 'bed_files', basename => "${sample}_peaks.narrowPeak", type => 'txt', metadata => { $sample_key => $sample });
                my $broadPeak_file  = $self->output_file(output_key => 'bed_files', basename => "${sample}_peaks.broadPeak",  type => 'txt', metadata => { $sample_key => $sample });
                my $gappedPeak_file = $self->output_file(output_key => 'bed_files', basename => "${sample}_peaks.gappedPeak", type => 'txt', metadata => { $sample_key => $sample });
                my $R_script        = $self->output_file(output_key => 'r_file',    basename => "${sample}_model.r",          type => 'txt', metadata => { $sample_key => $sample });
                my $pdf_file        = $self->output_file(output_key => 'plot_file', basename => "${sample}_model.pdf",        type => 'any', metadata => { $sample_key => $sample });
                my @outfiles = ($xls_file, $narrowPeak_file, $broadPeak_file, $gappedPeak_file, $R_script, $pdf_file);
                
                my $macs2_options = "";
                foreach my $marker (keys %marker_options) {
                    if ($sample =~ /$marker/i) {
                        $macs2_options = $marker_options{$marker};
                    }
                }
                my $peak_file;
                if ($macs2_options =~ /broad/i) {
                    $peak_file = $broadPeak_file->path;
                }
                else {
                    $peak_file = $narrowPeak_file->path;
                }
                
                #output files when $macs2_options=~/-B|-bdg/:
                my $bedGraph1_file = $self->output_file(output_key => 'bdg_files', basename => "${sample}_treat_pileup.bdg",   type => 'txt', metadata => { $sample_key => $sample });
                my $bedGraph2_file = $self->output_file(output_key => 'bdg_files', basename => "${sample}_control_lambda.bdg", type => 'txt', metadata => { $sample_key => $sample });
                #my $bedGraph3_file = $self->output_file(output_key => 'bdg_files', basename => "${sample}_treat_pvalue.bdg", type => 'txt', metadata => {$sample_key => $sample});
                #my $bedGraph4_file = $self->output_file(output_key => 'bdg_files', basename => "${sample}_treat_qvalue.bdg", type => 'txt', metadata => {$sample_key => $sample});
                push(@outfiles, ($bedGraph1_file, $bedGraph2_file)); #,$bedGraph3_file,$bedGraph4_file));
                my $req      = $self->new_requirements(memory => 2000, time => 1);
                my $bam_path = $input_bam->path;
                my $this_cmd = "use VRPipe::Steps::macs_callpeak; VRPipe::Steps::macs_callpeak->run_macs2( macs2 => q[$macs2_exe], samtools => q[$samtools_exe], bedtools => q[$bedtools_exe], control => q[$control_bam], bamfile => q[$bam_path], macs2_options => q[$macs2_options], sample => q[$sample], peakfile => q[$peak_file], r_cmd_prefix => q[$r_cmd_prefix]);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@outfiles });
            }
        };
    }
    
    method outputs_definition {
        return {
            xls_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'called peaks in xls',
            ),
            bed_files => VRPipe::StepIODefinition->create(
                type            => 'txt',
                max_files       => -1,
                description     => 'peak locations and peak summits in bed formats',
                check_existence => 0,
            ),
            r_file => VRPipe::StepIODefinition->create(
                type            => 'txt',
                max_files       => -1,
                description     => 'R script produced by macs2',
                check_existence => 0,
            ),
            plot_file => VRPipe::StepIODefinition->create(
                type            => 'any',
                max_files       => -1,
                description     => 'pdf plot',
                check_existence => 0,
            ),
            bdg_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'bedGraph files',
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs ChIP-seq peak-calling using macs2 callpeak";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method run_macs2 (ClassName|Object $self: Str :$macs2!, Str :$samtools!, Str :$bedtools!, Str :$control!, Str :$bamfile!, Str :$macs2_options!, Str :$sample!, Str :$peakfile!, Str :$r_cmd_prefix!) {
        my $input_bam = VRPipe::File->get(path => file($bamfile));
        my $peak_file = VRPipe::File->get(path => file($peakfile));
        my $ctrl_str  = "";
        if ($control) {
            $ctrl_str = " -c $control";
        }
        my $cmd_line = "$macs2 callpeak$ctrl_str -t $bamfile -n $sample -B $macs2_options";
        $input_bam->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $peak_file->update_stats_from_disc(retries => 3);
        
        $input_bam->add_metadata({ num_of_peaks => $peak_file->num_records });
        
        my $R_script = $peak_file->dir . "/${sample}_model.r";
        if (-e "$R_script") {
            my $this_cmd = $r_cmd_prefix . " $R_script";
            $input_bam->disconnect;
            system($this_cmd) && $self->throw("failed to run [$this_cmd]");
        }
        
        #Proportion of fragments in peaks (FrIP): first-cut metric for the success of the immunoprecipitation
        my $samtools_filter;
        if ($input_bam->metadata->{paired}) {
            $samtools_filter = "-f 3";
        }
        else {
            $samtools_filter = "-F 4";
        }
        my $FRIP = `$samtools view -b $samtools_filter $bamfile | $bedtools intersect -abam stdin -b $peakfile -bed -wa -u -sorted | cut -f 4 | sed 's/#.*//' | sort | uniq | wc -l`;
        chomp($FRIP);
        $input_bam->add_metadata({ FRIP => $FRIP });
    
    }

}

1;
