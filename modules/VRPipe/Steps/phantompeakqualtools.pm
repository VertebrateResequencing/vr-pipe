
=head1 NAME

VRPipe::Steps::phantompeakqualtools - a step

=head1 DESCRIPTION

The step computes highly informative enrichment and quality measures for 
ChIP-seq, DNase-seq, FAIRE-seq, MNase-seq using cross-correlation analysis.  It
also produces three quality metrics: NSC, RSC, and PBC which are measures  of
signal-to-noise and library complexity. See: doi: 10.1534/g3.113.008680. The
estimated quality measures are then added to metadata of input bam files.

It is very important that input bams are pre-filtered for unmapped reads,  low
quality reads, multimapping and duplicate reads: 1) Large number of
multimapping reads can severly affect the phantom peak  coefficient. 2) The
phantom peak coefficient can be artificially good for data with high  PCR
bottlenecking.

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

class VRPipe::Steps::phantompeakqualtools extends VRPipe::Steps::r {
    around options_definition {
        return {
            %{ $self->$orig },
            phantompeak_script  => VRPipe::StepOption->create(description => 'The R script to compute the frag length, data quality characteristics based on cross-correlation analysis and/or peak calling', optional => 1, default_value => 'run_spp.R'),
            sample_metadata_key => VRPipe::StepOption->create(description => 'metadata key for sample name',                                                                                                  optional => 1, default_value => 'sample')
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
            my $sample_key         = $options->{sample_metadata_key};
            my $phantompeak_script = $options->{phantompeak_script};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => "$phantompeak_script",
                    version => '',
                    summary => "Rscript $phantompeak_script -c=\$input_bam -savp -out=\$outfile"
                )
            );
            
            foreach my $input_bam (@{ $self->inputs->{bam_files} }) {
                my $sample  = $input_bam->metadata->{$sample_key};
                my $profile = $self->output_file(output_key => 'text_files', basename => "$sample", type => 'txt', metadata => { $sample_key => $sample });
                my $logfile = $self->output_file(output_key => 'text_files', basename => "$sample.log", type => 'txt', metadata => { $sample_key => $sample });
                my $pdfile  = $self->output_file(output_key => 'plot_file', basename => "$sample.pdf", type => 'any', metadata => { $sample_key => $sample });
                
                my @outfiles = ($profile, $logfile, $pdfile);
                my $this_cmd = $self->r_cmd_prefix . " $phantompeak_script -rf -savp -c=" . $input_bam->path . " -odir=" . $profile->dir . " -out=$sample > $sample.log";
                my $req      = $self->new_requirements(memory => 2000, time => 1);
                $self->dispatch_wrapped_cmd('VRPipe::Steps::phantompeakqualtools', 'compute_and_add_metadata', [$this_cmd, $req, { output_files => \@outfiles }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            text_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'cross-correlation results',
            ),
            plot_file => VRPipe::StepIODefinition->create(
                type        => 'any',
                max_files   => -1,
                description => 'pdf cross-correlation plots',
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Computes ChIP-seq data quality measures based on relative phantom peak";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method compute_and_add_metadata (ClassName|Object $self: Str $cmd_line) {
        my ($in, $out_dir, $out) = $cmd_line =~ /-c=(\S+) -odir=(\S+) -out=(\S+)/;
        my $input_bam = VRPipe::File->get(path => $in);
        my $out_file  = VRPipe::File->get(path => "$out_dir/$out");
        
        $input_bam->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $fh   = $out_file->openr;
        my $line = <$fh>;
        chomp($line);
        my @scores = split(/\t/, $line);
        $self->throw("file $out_dir/$out does not contain 11 tab delimited columns") unless ($#scores == 10);
        $fh->close;
        
        if ($input_bam->metadata->{reads}) {
            my $bam_reads = $input_bam->metadata->{reads};
            $self->throw("phantompeakqualtools processed $scores[1] reads while input file $in contains $bam_reads reads. Make sure input bam does not include multiple-alignments and is filtered for -F 0x0204.") unless ($bam_reads == $scores[1]);
        } #This condition may be too stringent, so the step almost always is preceded by a bam filtering step.
        
        my $meta = {};
        my @keys = qw(Filename numReads estFragLen corr_estFragLen phantomPeak corr_phantomPeak argmin_corr min_corr NSC RSC QualityTag);
        for my $i (0 .. $#scores) {
            $meta->{"$keys[$i]"} = $scores[$i];
        }
        $input_bam->add_metadata($meta);
    }

}

1;
