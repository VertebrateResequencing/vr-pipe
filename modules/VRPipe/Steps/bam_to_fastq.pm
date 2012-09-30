
=head1 NAME

VRPipe::Steps::bam_to_fastq - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::bam_to_fastq with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { ignore_read_ordering => VRPipe::StepOption->create(description => 'Do not require the input BAM to be name sorted. Output paired fastq are not necessarily in the same order nor guaranteed to have equal number of reads.', optional => 1, default_value => 0) };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => '1 or more BAM files. If the ignore_read_ordering option not set (the default) then BAMs must be name sorted.',
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
            my $ignore_read_ordering = $options->{ignore_read_ordering};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            my $bam_id = 0;
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
                $basename = "${basename}_$bam_id";
                $bam_id++;
                
                my $out_spec;
                my @fastqs;
                if ($paired) {
                    my $fastq = $self->output_file(
                        output_key => 'fastq_files',
                        basename   => "$basename.1.fastq",
                        type       => 'fq',
                        metadata   => {
                            %$fastq_meta,
                            reads  => $meta->{forward_reads},
                            paired => 1
                        }
                    );
                    my $reverse = $self->output_file(
                        output_key => 'fastq_files',
                        basename   => "$basename.2.fastq",
                        type       => 'fq',
                        metadata   => {
                            %$fastq_meta,
                            reads  => $meta->{reverse_reads},
                            paired => 2
                        }
                    );
                    @fastqs = ($fastq, $reverse);
                    
                    $fastq->add_metadata({ mate => $reverse->path->stringify });
                    $reverse->add_metadata({ mate => $fastq->path->stringify });
                    
                    $out_spec = 'forward => q[' . $fastq->path . '], reverse => q[' . $reverse->path . ']';
                }
                else {
                    my $fastq = $self->output_file(
                        output_key => 'fastq_files',
                        basename   => "$basename.0.fastq",
                        type       => 'fq',
                        metadata   => {
                            %$fastq_meta,
                            reads           => $meta->{reads},
                            bases           => $meta->{bases},
                            avg_read_length => sprintf("%0.2f", $meta->{reads} / $meta->{bases}),
                            paired          => 0
                        }
                    );
                    @fastqs = ($fastq);
                    
                    $out_spec = 'single => q[' . $fastq->path . ']';
                }
                
                my $this_cmd = "use VRPipe::Steps::bam_to_fastq; VRPipe::Steps::bam_to_fastq->bam_to_fastq(bam => q[$source_bam], $out_spec, ignore_read_ordering => $ignore_read_ordering);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@fastqs });
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
    
    method bam_to_fastq (ClassName|Object $self: Str|File :$bam!, Str|File :$forward?, Str|File :$reverse?, Str|File :$single?, Bool :$ignore_read_ordering = 0) {
        if ((defined $forward ? 1 : 0) + (defined $reverse ? 1 : 0) == 1) {
            $self->throw("When forward is used, reverse is required, and vice versa");
        }
        if (!$forward && !$single) {
            $self->throw("At least one of single or forward+reverse are required");
        }
        
        my $in_file = VRPipe::File->get(path => $bam);
        my @out_files;
        my @out_fhs;
        my $i = -1;
        foreach my $fq_path ($forward, $reverse, $single) {
            $i++;
            next unless $fq_path;
            push(@out_files, VRPipe::File->get(path => $fq_path));
            $out_fhs[$i] = $out_files[-1]->openw;
        }
        
        # parse bam file, convert to fastq, using OQ values if present
        my $pars = VRPipe::Parser->create('bam', { file => $in_file });
        $in_file->disconnect;
        my $pr = $pars->parsed_record();
        $pars->get_fields('QNAME', 'FLAG', 'SEQ', 'QUAL', 'OQ');
        my %pair_data;
        my ($pair_count, $single_count, $total_reads_parsed, $forward_count, $reverse_count) = (0, 0, 0, 0, 0);
        while ($pars->next_record()) {
            my $qname = $pr->{QNAME};
            my $flag  = $pr->{FLAG};
            
            my $seq  = $pr->{SEQ};
            my $oq   = $pr->{OQ};
            my $qual = $oq eq '*' ? $pr->{QUAL} : $oq;
            if ($pars->is_reverse_strand($flag)) {
                $seq = reverse($seq);
                $seq =~ tr/ACGTacgt/TGCAtgca/;
                $qual = reverse($qual);
            }
            
            if ($forward && $pars->is_sequencing_paired($flag) && $ignore_read_ordering) {
                my $first = $pars->is_first($flag);
                my $fh = $first ? $out_fhs[0] : $out_fhs[1];
                if ($first) {
                    print $fh '@', $qname, "/1\n", $seq, "\n+\n", $qual, "\n";
                    $forward_count++;
                }
                else {
                    print $fh '@', $qname, "/2\n", $seq, "\n+\n", $qual, "\n";
                    $reverse_count++;
                }
            }
            elsif ($forward && $pars->is_sequencing_paired($flag)) {
                my $key = $pars->is_first($flag) ? 'forward' : 'reverse';
                $pair_data{$key} = [$qname, $seq, $qual];
                my $f = $pair_data{forward} || [''];
                my $r = $pair_data{reverse} || [''];
                if ($f->[0] eq $r->[0]) {
                    my $fh = $out_fhs[0];
                    print $fh '@', $f->[0], "/1\n", $f->[1], "\n+\n", $f->[2], "\n";
                    $fh = $out_fhs[1];
                    print $fh '@', $r->[0], "/2\n", $r->[1], "\n+\n", $r->[2], "\n";
                    $pair_count++;
                }
            }
            elsif ($single) {
                my $fh = $out_fhs[2];
                print $fh '@', $qname, "\n", $seq, "\n+\n", $qual, "\n";
                $single_count++;
            }
            
            $total_reads_parsed++;
        }
        
        foreach my $of (@out_files) {
            $of->close;
        }
        
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
            # Instead we just confirm that the number of reads we wanted to
            # write out could actually be read back again
            my $expected_reads = $fq_meta->{paired} ? $pair_count : $single_count;
            if ($fq_meta->{paired} && $ignore_read_ordering) {
                $expected_reads = $fq_meta->{paired} == 1 ? $forward_count : $reverse_count;
            }
            unless ($expected_reads == $these_reads) {
                $self->throw("Made fastq file " . $out_file->path . " but there were only $these_reads reads instead of $expected_reads");
            }
            $extra_meta{ $out_file->id }->{reads} = $these_reads;
            
            # we can't really do much checking for base counts, so will assume
            # that if the above is fine, bases must be fine
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
