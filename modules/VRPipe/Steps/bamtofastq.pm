
=head1 NAME

VRPipe::Steps::bamtofastq - a step

=head1 DESCRIPTION

Runs the biobambam's bamtofastq executable to extract sequences from a BAM 
file into fastq

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Steps::bamtofastq with VRPipe::StepRole {
    use VRPipe::Parser;
    use POSIX;
    
    method options_definition {
        return {
            bamtofastq_exe   => VRPipe::StepOption->create(description => 'path to bamtofastq executable',                                                                                                                                                      optional => 1, default_value => 'bamtofastq'),
            bamtofastq_opts  => VRPipe::StepOption->create(description => 'bamtofastq options (excluding arguments that set input/output file names)',                                                                                                          optional => 1, default_value => 'gz=1 exclude=SECONDARY,SUPPLEMENTARY,QCFAIL'),
            fastqcheck_exe   => VRPipe::StepOption->create(description => 'path to fastqcheck executable',                                                                                                                                                      optional => 1, default_value => 'fastqcheck'),
            fastq_chunk_size => VRPipe::StepOption->create(description => 'size of output fastq chunks in base pair. When a value other than zero is specified, bamtofastq split option is used to generate the output fastqs in chunks of the specified size', optional => 1, default_value => '1000000000'),
        };
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
                    avg_read_length  => 'the average length of reads',
                    library          => 'library name',
                    sample           => 'sample name',
                    center_name      => 'center name',
                    platform         => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study            => 'name of the study, put in the DS field of the RG header line',
                    optional         => ['lane', 'library', 'sample', 'center_name', 'platform', 'study', 'mean_insert_size']
                }
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self            = shift;
            my $options         = $self->options;
            my $bamtofastq_exe  = $options->{bamtofastq_exe};
            my $bamtofastq_opts = $options->{bamtofastq_opts};
            my $fastqcheck_exe  = $options->{fastqcheck_exe};
            my $chunk_size      = $options->{fastq_chunk_size};
            if ($bamtofastq_opts =~ /\s(F|F2|S|O|O2|filename|fasta|split)=/) {
                $self->throw("bamtofastq_opts should not include F=, F2=, S=, O=, O2=, filename=, fasta= or split= options!");
            }
            $self->throw("Negative values cannot be used for fastq_chunk_size!") unless $chunk_size >= 0;
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bamtofastq',
                    version => VRPipe::StepCmdSummary->determine_version($bamtofastq_exe . ' --version', '^This.+version (.+)\.$'),
                    summary => "bamtofastq $bamtofastq_opts filename=\$bam_file split=\$reads_num F=\$fastq_1 F2=\$fastq_2 S=\$fastq_M 2>\$log_file"
                )
            );
            
            my $fastq_suffix = '';
            if ($bamtofastq_opts =~ /gz=1/) {
                $fastq_suffix = ".gz";
            }
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $source_bam = $bam->path->stringify;
                my @fq_files = VRPipe::Steps::bamtofastq->split_fastq_outs(split_dir => $self->output_root, bam => $source_bam, chunk_size => $chunk_size, suffix => $fastq_suffix);
                my @outfiles;
                foreach my $out (@fq_files) {
                    push(@outfiles, $self->output_file(output_key => 'fastq_files', basename => $out->basename, type => 'fq'));
                }
                
                my $basename = $bam->basename;
                $basename =~ s/\.bam//;
                my $out_log = $self->output_file(temporary => 1, basename => "$basename.log", type => 'txt');
                push(@outfiles, $out_log);
                my @fq_paths = map { $_->path } @fq_files;
                
                my $this_cmd = "use VRPipe::Steps::bamtofastq; VRPipe::Steps::bamtofastq->bamtofastq(bam => q[$source_bam], fastqs => [qw(@fq_paths)], bamtofastq_exe => q[$bamtofastq_exe], bamtofastq_opts => q[$bamtofastq_opts], fastqcheck_exe => q[$fastqcheck_exe], chunk_size => q[$chunk_size]);";
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
                    optional         => ['mate', 'lane', 'library', 'sample', 'center_name', 'platform', 'study', 'insert_size', 'mean_insert_size']
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
    
    method bamtofastq (ClassName|Object $self: Str|File :$bam!, ArrayRef[Str|File] :$fastqs!, Str|File :$bamtofastq_exe, Str :$bamtofastq_opts, Str|File :$fastqcheck_exe, Int :$chunk_size!) {
        my $in_file  = VRPipe::File->get(path => $bam);
        my $bam_meta = $in_file->metadata;
        my $basename = $in_file->basename;
        $basename =~ s/\.bam//;
        
        #determine the number of reads per fastq split
        my $splits         = $self->_splits($bam, $chunk_size);
        my $read_lengths   = $bam_meta->{avg_read_length};
        my $seqs_per_split = floor($chunk_size / $read_lengths);
        my $split_args     = '';
        if ($splits != 1) {
            $split_args = " split=$seqs_per_split";
        }
        
        #currently bamtofastq appends .gz to filenames when split option is used,
        #so we only explicitly add the suffix when fastq files are not chunked:
        my $suffix = '';
        if ($bamtofastq_opts =~ /gz=1/ && $splits == 1) {
            $suffix = ".gz";
        }
        
        my $out_dir   = VRPipe::File->get(path => $fastqs->[0])->dir;
        my $logfile   = "$out_dir/$basename.log";
        my $out_param = "F=$out_dir/${basename}_1.fastq$suffix F2=$out_dir/${basename}_2.fastq$suffix S=$out_dir/${basename}_M.fastq$suffix";
        my $cmd       = "$bamtofastq_exe $bamtofastq_opts $out_param$split_args filename=$bam 2>$logfile";
        $in_file->disconnect;
        system($cmd) && $self->throw("failed to run [$cmd]");
        
        # Create additional fastq files if necessary, or get rid of the unneeded empty files
        #Note: bamtofastq may produce empty files (e.g. empty _M.fq for when paired=1) which
        #we remove if they are not defined in method split_fastq_outs. It may also produce
        #extra files e.g. non-empty _M.fq even if paired=1 when the bam contains both paired
        #and single-ended reads.
        my @out_files;
        my @fq_files = @$fastqs;
        while ($splits) {
            my @out_pairs = ();
            if ($bam_meta->{paired}) {
                my $forward = shift @fq_files;
                my $reverse = shift @fq_files;
                push(@out_pairs, VRPipe::File->get(path => $forward));
                push(@out_pairs, VRPipe::File->get(path => $reverse));
                my $path = $forward;
                $path =~ s/_1.fastq/_M.fastq/;
                add_or_rm_fastq($self, $path, $forward, \@out_pairs);
                push(@out_files, @out_pairs);
            }
            else {
                my $single = shift @fq_files;
                push(@out_pairs, VRPipe::File->get(path => $single));
                foreach my $n (1, 2) {
                    my $path = $single;
                    $path =~ s/_M.fastq/_$n.fastq/;
                    add_or_rm_fastq($self, $path, $single, \@out_pairs);
                    push(@out_files, @out_pairs);
                }
            }
            $splits--;
        }
        
        # check logfile read count agrees with bam meta
        my $log = VRPipe::File->get(path => "$logfile");
        $log->update_stats_from_disc;
        my $logh = $log->openr;
        my @lines;
        while (my $line = <$logh>) {
            unless ($line !~ /^\+$/) {
                $self->throw("$logfile suggests that bam file contains unmatched pair mates: $bam");
            }
            push(@lines, $line);
        }
        $logh->close;
        my $last = pop(@lines);
        unless ($last =~ /MemUsage/ && $last =~ /wall clock time/) { $self->throw("$logfile does not contain runtime stats") }
        $last = pop(@lines);
        my $total_reads_parsed = $1 if $last =~ /\[V\] (\S+)/;
        
        unless ($total_reads_parsed) {
            foreach my $out_file (@out_files) {
                $out_file->unlink;
            }
            $self->throw("Failed to get sequence count from $logfile, or count is zero");
        }
        unless ($total_reads_parsed == $bam_meta->{reads}) {
            foreach my $out_file (@out_files) {
                $out_file->unlink;
            }
            $self->throw("$logfile says we parsed $total_reads_parsed instead of bam meta $bam_meta->{reads} reads");
        }
        
        # check the fastq files are as expected
        my %extra_meta;
        foreach my $out_file (@out_files) {
            # get read and base counts from fastqcheck
            my ($these_reads, $these_bases) = (0, 0);
            my $path = $out_file->path;
            my $pipe;
            if ($bamtofastq_opts =~ /gz=1/) {
                $pipe = "gunzip -c $path | $fastqcheck_exe |";
            }
            else {
                $pipe = "$fastqcheck_exe $path |";
            }
            open(my $fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
            while (<$fh>) {
                # should be first line
                unless (/(\S+) sequences, (\S+) total length/) {
                    foreach my $out_file (@out_files) {
                        $out_file->unlink;
                    }
                    $self->throw("Could not get reads and bases from $pipe : $_");
                }
                ($these_reads, $these_bases) = ($1, $2);
                last;
            }
            close($fh);
            
            # for expected_reads we do not use the 'reads' metadata since that
            # is based on bamcheck forward/reverse reads, which is based on
            # flags, so if flags are wrong, or reads have been removed from
            # input bam, either way resulting in mismatch of forward and reverse
            # reads, we will have discarded some reads.
            $out_file->update_stats_from_disc;
            my $output_lines = $out_file->lines;
            unless ($output_lines == $these_reads) {
                foreach my $out_file (@out_files) {
                    $out_file->unlink;
                }
                $self->throw("$path record count $output_lines disagrees with fastqcheck count $these_reads");
            }
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
        
        return 1;
    }
    
    sub add_or_rm_fastq {
        # Create additional files if necessary, or get rid of the unneeded empty files
        my ($self, $path, $based_on, $out_files) = @_;
        
        #delete the extra_file=$path if it exists and is empty. keep it if it contains data
        if (-e $path) {
            my $pipe;
            if ($path =~ /\.gz$/) {
                $pipe = "(gunzip -c $path | wc -l; echo $path) | xargs -n2 |";
            }
            else {
                $pipe = "wc -l $path |";
            }
            
            my $file_size;
            open(my $fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
            while (<$fh>) {
                unless (/(\S+) $path/) {
                    foreach my $out_file (@$out_files) {
                        $out_file->unlink;
                    }
                    $self->throw("Couldn't get count from $pipe : $_");
                }
                $file_size = $1;
                last;
            }
            if ($file_size == 0) {
                unless (system("rm $path") == 0) {
                    foreach my $out_file (@$out_files) {
                        $out_file->unlink;
                    }
                    $self->throw("failed to rm $path");
                }
                return;
            }
        }
        else {
            return;
        }
        
        # figure out what stepstate we are for and add the fastq as an output file
        my @existing_stepoutputfiles = VRPipe::StepOutputFile->search({ file => VRPipe::File->get(path => $based_on)->id, output_key => 'fastq_files' });
        my %stepstates = map { $_->id => $_ } @existing_stepoutputfiles;
        if (keys %stepstates != 1) {
            foreach my $out_file (@$out_files) {
                $out_file->unlink;
            }
            $self->throw("Could not get unique stepstate for $based_on");
        }
        
        my $fastq_meta = VRPipe::File->get(path => $based_on)->metadata;
        if ($path =~ /_M.fastq/) {
            $fastq_meta->{paired} = 0;
            delete $fastq_meta->{mate};
        }
        elsif ($path =~ /_1.fastq/) {
            $fastq_meta->{paired} = 1;
            my $mate = $path;
            $mate =~ s/_1.fastq/_2.fastq/;
            $fastq_meta->{mate} = $mate;
        }
        else {
            $fastq_meta->{paired} = 2;
            my $mate = $path;
            $mate =~ s/_2.fastq/_1.fastq/;
            $fastq_meta->{mate} = $mate;
        }
        
        # create the extra_file=$path in dtabase and define it as output file
        my $extra_file = VRPipe::File->create(path => $path, metadata => {%$fastq_meta});
        
        my $step_state = $existing_stepoutputfiles[0]->stepstate;
        VRPipe::StepOutputFile->create(
            file       => $extra_file->id,
            stepstate  => $step_state,
            output_key => 'fastq_files',
        );
        push(@$out_files, $extra_file);
    }
    
    method _splits (ClassName|Object $self: Str|File $bam, Int $chunk_size) {
        my $bam_meta       = VRPipe::File->get(path => $bam)->metadata;
        my $read_lengths   = $bam_meta->{avg_read_length};
        my $seqs           = $bam_meta->{reads};
        my $seqs_per_split = floor($chunk_size / $read_lengths);
        
        if ($seqs_per_split != 0) {
            return ceil($seqs / $seqs_per_split / 2);
        }
        else {
            return 1;
        }
    }
    
    method split_fastq_outs (ClassName|Object $self: Str|Dir :$split_dir!, Str|File :$bam!, Int :$chunk_size!, Str :$suffix!) {
        my $splits = $self->_splits($bam, $chunk_size);
        
        my $bam_file   = VRPipe::File->get(path => $bam);
        my $meta       = $bam_file->metadata;
        my $paired     = $meta->{paired};
        my $source_bam = $bam_file->path->stringify;
        my $fastq_meta = { source_bam => $source_bam };
        foreach my $key (qw(lane insert_size mean_insert_size library sample center_name platform study split_sequence population continent chrom from to individual project species public_name sample_accession_number)) {
            if (defined $meta->{$key}) {
                $fastq_meta->{$key} = $meta->{$key};
            }
        }
        my @outs;
        my $basename = $bam_file->basename;
        $basename =~ s/\.bam//;
        
        #We expect to generate _1.fq & _2.fq for paired reads and _M.fq for single-ended reads.
        #if more files are produced by bamtofastq we define them later in method add_or_rm_fastq
        if ($splits == 1) {
            if ($paired) {
                my $forward = VRPipe::File->create(
                    path => file($split_dir, "${basename}_1.fastq$suffix"),
                    type => 'fq',
                    metadata => { %$fastq_meta, paired => 1, chunk => 1 }
                );
                my $reverse = VRPipe::File->create(
                    path => file($split_dir, "${basename}_2.fastq$suffix"),
                    type => 'fq',
                    metadata => { %$fastq_meta, paired => 2, chunk => 1 }
                );
                @outs = ($forward, $reverse);
                $forward->add_metadata({ mate => $reverse->path->stringify });
                $reverse->add_metadata({ mate => $forward->path->stringify });
            }
            else {
                my $single = VRPipe::File->create(
                    path => file($split_dir, "${basename}_M.fastq$suffix"),
                    type => 'fq',
                    metadata => { %$fastq_meta, paired => 0, chunk => 1 }
                );
                @outs = ($single);
            }
        }
        else {
            for my $split_num ("000000" .. sprintf("%06d", $splits - 1)) {
                if ($paired) {
                    my $forward = VRPipe::File->create(
                        path => file($split_dir, "${basename}_1.fastq_$split_num$suffix"),
                        type => 'fq',
                        metadata => { %$fastq_meta, paired => 1, chunk => $split_num }
                    );
                    my $reverse = VRPipe::File->create(
                        path => file($split_dir, "${basename}_2.fastq_$split_num$suffix"),
                        type => 'fq',
                        metadata => { %$fastq_meta, paired => 2, chunk => $split_num }
                    );
                    push(@outs, $forward, $reverse);
                    $forward->add_metadata({ mate => $reverse->path->stringify });
                    $reverse->add_metadata({ mate => $forward->path->stringify });
                }
                else {
                    my $single = VRPipe::File->create(
                        path => file($split_dir, "${basename}_M.fastq_$split_num$suffix"),
                        type => 'fq',
                        metadata => { %$fastq_meta, paired => 0, chunk => $split_num }
                    );
                    push(@outs, $single);
                }
            }
        }
        
        return @outs;
    }
}

1;
