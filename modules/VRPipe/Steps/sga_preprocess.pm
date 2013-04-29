
=head1 NAME

VRPipe::Steps::sga_preprocess - a step

=head1 DESCRIPTION

Prepare nated fastq files for assembly with sga Requires that the fastq files
are suffixed with _[0|1|M].fastq(.gz) as per bam2fastq executable output, where
1=forward reads, 2=reverse, M=single (Currently we reject any single reads from
the process)

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

# Usage: sga preprocess [OPTION] --quality-scale=STR READS1 READS2 ...
# Prepare READS1, READS2, ... data files for assembly
# If pe-mode is turned on (pe-mode=1) then if a read is discarded its pair will be discarded as well.
#
#       --help                           display this help and exit
#       -v, --verbose                    display verbose output
#       -o, --out=FILE                   write the reads to FILE (default: stdout)
#           --phred64                    the input reads are phred64 scaled. They will be converted to phred33.
#       -p, --pe-mode=INT                0 - do not treat reads as paired (default)
#                                        1 - reads are paired with the first read in READS1 and the second
#                                        read in READS2. The paired reads will be interleaved in the output file
#                                        2 - reads are paired and the records are interleaved within a single file.
#       -q, --quality-trim=INT           perform Heng Li's BWA quality trim algorithm.
#                                        Reads are trimmed according to the formula:
#                                        argmax_x{\sum_{i=x+1}^l(INT-q_i)} if q_l<INT
#                                        where l is the original read length.
#       -f, --quality-filter=INT         discard the read if it contains more than INT low-quality bases.
#                                        Bases with phred score <= 3 are considered low quality. Default: no filtering.
#                                        The filtering is applied after trimming so bases removed are not counted.
#       -r, --remove-adapter-fwd=STRING
#       -c, --remove-adapter-rev=STRING  Remove the adapter STRING from input reads.
#                                        Do not use this option if you are planning to use the BCR algorithm for indexing.
#       -m, --min-length=INT             discard sequences that are shorter than INT
#                                        this is most useful when used in conjunction with --quality-trim. Default: 40
#       -h, --hard-clip=INT              clip all reads to be length INT. In most cases it is better to use
#                                        the soft clip (quality-trim) option.
#       --permute-ambiguous              Randomly change ambiguous base calls to one of possible bases.
#                                        For example M will be changed to A or C. If this option is not specified, the
#                                        entire read will be discarded.
#       -s, --sample=FLOAT               Randomly sample reads or pairs with acceptance probability FLOAT.
#       --dust                           Perform dust-style filtering of low complexity reads. If you are performing
#                                        de novo genome assembly, you probably do not want this.
#       --dust-threshold=FLOAT           filter out reads that have a dust score higher than FLOAT (default: 4.0).
#                                        This option implies --dust
#       --suffix=SUFFIX                  append SUFFIX to each read ID

use VRPipe::Base;

class VRPipe::Steps::sga_preprocess with VRPipe::StepRole {
    method options_definition {
        return {
            sga_preprocess_options        => VRPipe::StepOption->create(description => 'options to sga preprocess',                             optional => 1, default_value => '--min-length=75'),
            sga_exe                       => VRPipe::StepOption->create(description => 'path to your sga executable',                           optional => 1, default_value => 'sga'),
            sga_preprocess_compress_fastq => VRPipe::StepOption->create(description => 'compress the fastq output of sga preprocess (boolean)', optional => 1, default_value => 1)
        };
    }
    
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->create(type => 'fq', max_files => -1, description => 'fastq files to be indexed') };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $sga_exe  = $options->{sga_exe};
            my $sga_opts = $options->{sga_preprocess_options};
            my $compress = $options->{sga_preprocess_compress_fastq};
            if ($sga_opts =~ /preprocess/) {
                $self->throw("sga_preprocess_options should not include the preprocess subcommand");
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'sga', version => VRPipe::StepCmdSummary->determine_version($sga_exe, '^Version: (.+)$'), summary => 'sga preprocess ' . $sga_opts . ' $fastq_file(s)'));
            
            my %fastqs;
            foreach my $fq (@{ $self->inputs->{fastq_files} }) {
                my $meta = $fq->metadata;
                next unless $meta->{paired};
                my $basename = $fq->basename;
                $basename =~ s/_(1|2|M)\.(fq|fastq)(\.gz)?$/\.processed.fq/;
                if ($compress) {
                    $basename .= '.gz';
                }
                if ($meta->{paired} == 1) {
                    unshift @{ $fastqs{$basename} }, $fq;
                }
                else {
                    push @{ $fastqs{$basename} }, $fq;
                }
            }
            
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $fq (keys %fastqs) {
                my @fqs          = map { $_->path } @{ $fastqs{$fq} };
                my $meta         = $self->common_metadata($fastqs{$fq});
                my $processed_fq = $self->output_file(output_key => 'preprocessed_fastq_files', basename => $fq, type => 'fq', metadata => $meta);
                my $write_to     = $compress ? '| gzip -c >' : '>';
                my $cmd          = qq[$sga_exe preprocess $sga_opts ] . join(' ', @fqs) . " $write_to " . $processed_fq->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::sga_preprocess', 'preprocess_and_check', [$cmd, $req, { output_files => [$processed_fq] }]);
            }
        };
    }
    
    method outputs_definition {
        return { preprocessed_fastq_files => VRPipe::StepIODefinition->create(type => 'fq', description => 'the preprocessed fastq files', max_files => -1) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Prepare fastq files for assembly with sga";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method preprocess_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($out_path) = $cmd_line =~ /> (\S+)$/;
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $out_fq = VRPipe::File->get(path => $out_path);
        $out_fq->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        $out_fq->update_stats_from_disc(retries => 3);
        
        my $reads = $out_fq->num_records;
        if ($reads > 0) {
            $out_fq->add_metadata({ reads => $reads });
            return 1;
        }
        else {
            $out_fq->unlink;
            $self->throw("cmd [$cmd_line] failed because no reads were generated in the output fastq file - nothing passes you preprocessing thresholds??");
        }
    }
}

1;
