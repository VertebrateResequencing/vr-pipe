
=head1 NAME

VRPipe::Steps::run_bam2fastq - a step

=head1 DESCRIPTION

Runs the bam2fastq executable to extract sequences from a BAM file into fastq

=head1 AUTHOR

Chris Joyce <ch5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::run_bam2fastq with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { 
            bam2fastq_exe => VRPipe::StepOption->create(description => 'path to bam2fastq executable', optional => 1, default_value => 'bam2fastq'),
            bam2fastq_opts => VRPipe::StepOption->create(description => 'bam2fastq options excluding --o', optional => 1),
        }
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => '1 or more BAM files',
                metadata    => {
                    lane             => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    bases            => 'total number of base pairs',
                    reads            => 'total number of reads (sequences)',
                    forward_reads    => 'number of forward reads',
                    reverse_reads    => 'number of reverse reads',
                    paired           => '0=single ended reads only; 1=paired end reads present',
                    mean_insert_size => 'mean insert size (0 if unpaired)',
                    library          => 'library name',
                    sample           => 'sample name',
                    center_name      => 'center name',
                    platform         => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study            => 'name of the study, put in the DS field of the RG header line',
                    optional         => ['library', 'sample', 'center_name', 'platform', 'study', 'mean_insert_size']
                }
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self                 = shift;
            my $options              = $self->options;
            my $bam2fastq_exe = $options->{bam2fastq_exe};
            my $bam2fastq_opts = $options->{bam2fastq_opts};
            $bam2fastq_opts  .= ' -f ';   # force overwrite
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $meta = $bam->metadata;
                next unless $meta->{reads};
                my $paired = $meta->{paired};
                
                my $source_bam = $bam->path->stringify;
                my $fastq_meta = { source_bam => $source_bam };
                foreach my $key (qw(lane insert_size mean_insert_size library sample center_name platform study split_sequence population continent chrom from to)) {
                    if (defined $meta->{$key}) {
                        $fastq_meta->{$key} = $meta->{$key};
                    }
                }
                my $basename = $bam->basename;
                $basename =~ s/\.bam//;
                
                my $out_spec;
                my @outfiles;
                if ($paired) {
                    my $fastq = $self->output_file(
                        output_key => 'fastq_files',
                        basename   => "${basename}_1.fastq",
                        type       => 'fq',
                        metadata   => {
                            %$fastq_meta,
                            reads  => $meta->{forward_reads},
                            paired => 1
                        }
                    );
                    my $reverse = $self->output_file(
                        output_key => 'fastq_files',
                        basename   => "${basename}_2.fastq",
                        type       => 'fq',
                        metadata   => {
                            %$fastq_meta,
                            reads  => $meta->{reverse_reads},
                            paired => 2
                        }
                    );
                    @outfiles = ($fastq, $reverse);
                    
                    $fastq->add_metadata({ mate => $reverse->path->stringify });
                    $reverse->add_metadata({ mate => $fastq->path->stringify });
                    
                    $out_spec = 'forward => q[' . $fastq->path . '], reverse => q[' . $reverse->path . ']';
                }
                else {
                    my $fastq = $self->output_file(
                        output_key => 'fastq_files',
                        basename   => "${basename}_M.fastq",
                        type       => 'fq',
                        metadata   => {
                            %$fastq_meta,
                            reads           => $meta->{reads},
                            bases           => $meta->{bases},
                            avg_read_length => sprintf("%0.2f", $meta->{reads} / $meta->{bases}),
                            paired          => 0
                        }
                    );
                    @outfiles = ($fastq);
                    
                    $out_spec = 'single => q[' . $fastq->path . ']';
                }
                
                my $out_log = $self->output_file(output_key => 'bam2fastq_logs', basename => "$basename.log", type => 'txt');
                push (@outfiles, $out_log);

                my $this_cmd = "use VRPipe::Steps::run_bam2fastq; VRPipe::Steps::run_bam2fastq->run_bam2fastq(bam => q[$source_bam], $out_spec, bam2fastq_exe => q[$bam2fastq_exe], bam2fastq_opts => q[$bam2fastq_opts]);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@outfiles });
            }
        };
    }
    
    method outputs_definition {
        return {
            fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => -1,
                description => '1 or more fastq files',
                metadata    => {
                    lane             => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    bases            => 'total number of base pairs',
                    reads            => 'total number of reads (sequences)',
                    avg_read_length  => 'the average length of reads',
                    paired           => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                    mate             => 'if paired, the path to the fastq that is our mate',
                    source_bam       => 'path of the bam file this fastq file was made from',
                    insert_size      => 'expected library insert size (0 if unpaired)',
                    mean_insert_size => 'calculated mean insert size (0 if unpaired)',
                    library          => 'library name',
                    sample           => 'sample name',
                    center_name      => 'center name',
                    platform         => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study            => 'name of the study',
                    optional         => ['mate', 'library', 'sample', 'center_name', 'platform', 'study', 'insert_size', 'mean_insert_size']
                }
            ),
            bam2fastq_logs => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'bam2fastq log',
                max_files   => -1
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Converts bam files to fastq files";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method run_bam2fastq (ClassName|Object $self: Str|File :$bam!, Str|File :$forward?, Str|File :$reverse?, Str|File :$single?, Str|File :$bam2fastq_exe, Str :$bam2fastq_opts) {

        if ((defined $forward ? 1 : 0) + (defined $reverse ? 1 : 0) == 1) {
            $self->throw("When forward is used, reverse is required, and vice versa");
        }
        if (!$forward && !$single) {
            $self->throw("At least one of single or forward+reverse are required");
        }
        
        my $in_file = VRPipe::File->get(path => $bam);
        my @out_files;
        my $out_dir;
        foreach my $fq_path ($forward, $reverse, $single) {
            next unless $fq_path;
            my $fq = VRPipe::File->get(path => $fq_path);
            push(@out_files, $fq);
            $out_dir = $fq->dir unless $out_dir; 
        }
        
        my $basename = $in_file->basename;
        $basename =~ s/\.bam//;
        my $logfile = "$out_dir/$basename.log";
        # '#' is a special char, (replaced with _1 and _2 for PE reads, _M for SE reads).
        my $outfile = "$out_dir/${basename}#.fastq";

        my $cmd = "$bam2fastq_exe $bam2fastq_opts -o $outfile $bam 2>$logfile";
        system($cmd) && $self->throw("failed to run [$cmd]");

        my $total_reads_parsed;
        my $log = VRPipe::File->get(path => "$logfile");
        $log->update_stats_from_disc;
        my $logh = $log->openr;
        while (my $line = <$logh>) {
            $total_reads_parsed = $1 if $line =~ /(\S+) sequences in the BAM file/;
        }
        $self->throw("Failed to get sequence count from log, or count is zero") unless $total_reads_parsed;
        
        my $bam_meta = $in_file->metadata;
        unless ($total_reads_parsed == $bam_meta->{reads}) {
            foreach my $out_file (@out_files) {
                $out_file->unlink;
            }
            $self->throw("When parsing the bam file we parsed $total_reads_parsed instead of $bam_meta->{reads} reads");
        }

        # check the fastq files are as expected
        my %extra_meta;
        foreach my $out_file (@out_files) {

            next if $out_file->basename =~ /\.log$/;
            $out_file->update_stats_from_disc(retries => 3);
            
            # check that sequence lengths match quality string lengths, and
            # get read and base counts
            my ($these_reads, $these_bases) = (0, 0);
            my $pars = VRPipe::Parser->create('fastq', { file => $out_file });
            $in_file->disconnect;
            my $pr = $pars->parsed_record;
            while ($pars->next_record()) {
                my $id       = $pr->[0];
                my $seq_len  = length($pr->[1]);
                my $qual_len = length($pr->[2]);
                unless ($seq_len == $qual_len) {
                    $out_file->unlink;
                    $self->throw("Made fastq file " . $out_file->path . " but sequence $id had mismatching sequence and quality lengths ($seq_len vs $qual_len)");
                }
                $these_reads++;
                $these_bases += $seq_len;
            }
            
            my $fq_meta = $out_file->metadata;

            # for expected_reads we do not use the 'reads' metadata since that
            # is based on bamcheck forward/reverse reads, which is based on
            # flags, so if flags are wrong, or reads have been removed from
            # input bam, either way resulting in mismatch of forward and reverse
            # reads, we will have discarded some reads.

            $extra_meta{ $out_file->id }->{reads} = $these_reads;
            $extra_meta{ $out_file->id }->{bases} = $these_bases;
        }
        
        # add/correct metadata
        foreach my $out_file (@out_files) {
            my $extra = $extra_meta{ $out_file->id } || next;
            $out_file->add_metadata({
                    bases           => $extra->{bases},
                    reads           => $extra->{reads},
                    avg_read_length => sprintf("%0.2f", $extra->{bases} / $extra->{reads})
                },
                replace_data => 1
            );
        }
    }
}

1;
