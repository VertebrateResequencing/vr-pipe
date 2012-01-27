use VRPipe::Base;

class VRPipe::Steps::qc_stats_exome with VRPipe::StepRole {
    use VRPipe::Parser;
    use VRPipe::Parser::fastq;
    use VRPipe::Utils::Math;
	use List::Util qw(min max sum);
	use Data::Dumper;
	
	method options_definition {
        return { load_intervals_dump_file => VRPipe::StepOption->get(description => 'intervals dump file for exome data', optional => 1),
                 bait_interval => VRPipe::StepOption->get(description => 'filename of baits intervals file, of the form 1 bait per line', optional => 1),
                 target_interval => VRPipe::StepOption->get(description => 'filename of targets intervals file, same format as bait_interval', optional => 1),
                 ref_fa => VRPipe::StepOption->get(description => 'filename of reference fasta file', optional => 1),
                 ref_fai  => VRPipe::StepOption->get(description => 'fai file of reference fasta file', optional => 1),
                 name_check => VRPipe::StepOption->get(description => 'boolean which if true means that the sanity check of the reference sequence names is performed', optional => 1, default_value => 1)
               };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', description => 'bam files', max_files => -1) };
    }
    
    method body_sub {
        return sub { 
            my $self = shift;
            my $options = $self->options;
            my $extra_args = '';
            foreach my $extra (qw(load_intervals_dump_file name_check)) {
                if ($options->{$extra}) {
                    $extra_args .= ", $extra => q[$options->{$extra}]";
                }
            }
            my $req = $self->new_requirements(memory => 25000, time => 15, cpus => 4);
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
                my $sample;
                my $meta = $bam->metadata;
                if ($meta->{sample}) {
                	$sample = $meta->{sample};
                }
                else {
             	    my @samtools_sample = split("\t", `samtools view -H $bam_path | grep SM`);
             	    for my $sam_out ( @samtools_sample ) {
             	    	if ( $sam_out =~ s/^SM:// ) {
             	    		$sample = $sam_out;
             	    	}
             	    }
                }
                my $stats = $self->output_file(output_key => 'stats_dump_files',
                                               basename => $bam->basename.'.stats.dump',
                                               type => 'txt',
                                               metadata => {sample => $sample});
                my $stats_path = $stats->path;
                my $cmd = "use VRPipe::Steps::qc_stats_exome; VRPipe::Steps::qc_stats_exome->bam_exome_qc_stats(bam => q[$bam_path], dump => q[$stats_path]$extra_args);";
                $self->dispatch_vrpipecode($cmd, $req, {output_files => [$stats]});
            }
        };
    }
    
    method outputs_definition {
        return { stats_dump_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'dump file to store result of bam_exome_qc_stats subroutine per bam file', max_files => -1, metadata => {sample => 'sample name'})};
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Generates statistics dump of exome bam file using either a pre-calculated intervals dump file or a reference fasta and target/bait intervals files";
    }
    method max_simultaneous {
        return 0;
    }
    
    method bam_statistics (ClassName|Object $self: Str|File $bam, HashRef $stats, HashRef $intervals) {
        my $bam_file = VRPipe::File->get(path => $bam);
        my $pars = VRPipe::Parser->create('bam', {file => $bam_file});
        $pars->get_fields('FLAG', 'RNAME', 'POS', 'CIGAR', 'ISIZE', 'SEQ', 'QUAL', 'SEQ_LENGTH', 'MAPPED_SEQ_LENGTH', 'NM');
        $bam_file->disconnect;
        my $parsed_record = $pars->parsed_record();
        my $first_sam_record = 1;
        my $max_insert = 1000;
        my $low_cvg = 2;
        my $current_rname;
        my %pileup;
		while ($pars->next_record()) {
            my $flag = $parsed_record->{FLAG};
            my $rname = $parsed_record->{RNAME};
            my $pos = $parsed_record->{POS};
            my $cigar = $parsed_record->{CIGAR};
            my $isize = $parsed_record->{ISIZE};
            my $seq = $parsed_record->{SEQ};
            my $qual = $parsed_record->{QUAL};
            my $seq_length = $parsed_record->{SEQ_LENGTH};
            my $mapped_seq_length = $parsed_record->{MAPPED_SEQ_LENGTH};
            my $num_mismatches = $parsed_record->{NM};
            my $gc_array = $$stats{gc_hist_unmapped_1};
            my $qual_array = $$stats{qual_scores_up};
            
            if ($first_sam_record) {
                $$stats{readlen} = length $qual;
                $$stats{paired} = $pars->is_sequencing_paired($flag);
                $$stats{gc_hist_unmapped_1} = [(0) x ($$stats{readlen} + 1)];
                $$stats{gc_hist_unmapped_2} = [(0) x ($$stats{readlen} + 1)];
                $$stats{gc_hist_mapped_1} = [(0) x ($$stats{readlen} + 1)];
                $$stats{gc_hist_mapped_2} = [(0) x ($$stats{readlen} + 1)];
                foreach my $a qw(gc_vs_bait_cvg gc_vs_target_cvg) {
                    $$stats{$a} = ();
                    foreach (0 .. 100) {
                        push @{$$stats{$a}}, {};
                    }
                }

                foreach my $a qw(qual_scores_1 qual_scores_2 qual_scores_up) {
                    $$stats{$a} = ();
                    foreach (1 .. $$stats{readlen}) {
                        push @{$$stats{$a}}, {};
                    }
                }

                $first_sam_record = 0;
            }
            $$stats{raw_reads}++; 
            $$stats{raw_bases} += length $qual;
            $$stats{reads_paired}++ if $pars->is_sequencing_paired($flag);
            $$stats{num_mismatches} += $num_mismatches unless $num_mismatches eq '*';
            
            if ($pars->is_mapped($flag)) {
                $$stats{reads_mapped}++;
                $$stats{bases_mapped} += $mapped_seq_length;
                $$stats{mapped_as_pair}++ if $pars->is_mapped_paired($flag);
                $$stats{insert_size_hist}{$isize}++ if ($pars->is_sequencing_paired($flag) and $isize > 0 and $isize < $max_insert);
                $$stats{bases_mapped_reads}+= $seq_length;
            
                unless ($pars->is_duplicate($flag)) {
                    $$stats{rmdup_reads_mapped}++;
                    $$stats{rmdup_bases_mapped} += $mapped_seq_length;
                }
            
                if ($pars->is_first($flag)) {
                    $gc_array = $$stats{gc_hist_mapped_1};
                }
                else {
                    $gc_array = $$stats{gc_hist_mapped_2}; 
                }
            }
            elsif ($pars->is_second($flag)) {
                $gc_array = $$stats{gc_hist_unmapped_2};
            }
            
            $gc_array->[$seq =~ tr/cgCG/cgCG/]++;

            $qual_array = $$stats{qual_scores_1} if $pars->is_first($flag);
            $qual_array = $$stats{qual_scores_2} if $pars->is_second($flag);
            my @qual_ints = VRPipe::Parser::fastq->qual_to_ints($qual); 
            @qual_ints = reverse @qual_ints if $pars->is_reverse_strand($flag);
            foreach my $i (0 .. $#qual_ints){
                $qual_array->[$i]{$qual_ints[$i]}++;
            }
            
            next unless $pars->is_mapped($flag);

            unless (exists $intervals->{target}{$rname} or exists $intervals->{bait}{$rname}) {
                $$stats{off_bait_bases} += $mapped_seq_length; 
                $$stats{off_target_bases} += $mapped_seq_length;
                next;
            }
            
            if ((!defined $current_rname) or $current_rname ne $rname) {
                if (defined $current_rname) {
                    $self->_update_bam_exome_qc_stats($intervals, \%pileup, $stats, $current_rname, $low_cvg);
                }

                $current_rname = $rname;
                # initialize the pileup hashes for the new reference sequence, which
                # means setting the value of each position to be zero
                foreach (keys %{$intervals}) {
                    $pileup{$_} = {};
                }  
                foreach my $interval_type ('bait', 'target') {
                    foreach my $interval_array ($intervals->{$interval_type}{$current_rname}) {
                        foreach my $interval (@{$interval_array}) {
                            foreach my $i ($interval->[0] .. $interval->[1]) {
                                $pileup{$interval_type}{$i} = 0;
                            }
                        }
                    }
                }

                # That's baits/targets initialized.  Now need to do the near bases.
                # We don't want a 'near base' to be actually in a target/bait, since doing so
                # would mean we're recording the same information twice, using extra memory
                foreach my $interval_type ('bait', 'target') {
                    foreach my $interval_array ($intervals->{$interval_type . '_near'}{$current_rname}) {
                        foreach my $interval (@{$interval_array}) {
                            foreach my $i ($interval->[0] .. $interval->[1]) {
                                $pileup{$interval_type . '_near'}{$i} = 0 unless exists $pileup{$interval_type}{$i};
                            }
                        }
                    }
                }
            }

            # only gather bait/target stats for reads which are not duplicates
            next if $pars->is_duplicate($flag);

            # update the number of reads hitting targets/baits and the target/bait coverage
            my $ref_pos = $pos;  
            my @numbers = split /[A-Z]/, $cigar;
            my @letters = split /[0-9]+/, $cigar;
            shift @letters;
            my %read_hits;
        
            foreach my $letter_index (0 ..  $#letters) {
                if ($letters[$letter_index] eq 'M') {
			        foreach my $i ($ref_pos .. $ref_pos + $numbers[$letter_index] - 1) {
                        foreach my $key ('bait', 'target') {
                            if (exists $pileup{$key}{$i}) {
                                $pileup{$key}{$i}++;
                                $read_hits{$key} = 1;
                                $read_hits{$key . '_near'} = 1;
                            }
                            elsif (exists $pileup{$key . '_near'}{$i}) {
                                $pileup{$key . '_near'}{$i}++;
                                $read_hits{$key . '_near'} = 1;
                            }
                        }

                        $$stats{off_bait_bases}++ unless (exists $pileup{bait}{$i} or exists $pileup{bait_near}{$i});
                        $$stats{off_target_bases}++ unless (exists $pileup{target}{$i} or exists $pileup{target_near}{$i});
                    } 

                    $ref_pos += $numbers[$letter_index];
                 }
                 elsif ($letters[$letter_index] eq 'D' or $letters[$letter_index] eq 'N') {
                     $ref_pos += $numbers[$letter_index];
                 }
                elsif ($letters[$letter_index] eq 'S' or $letters[$letter_index] eq 'H') {
                    $$stats{clip_bases} += $numbers[$letter_index];
                }
            }

            delete $read_hits{bait_near} if $read_hits{bait};
            delete $read_hits{target_near} if $read_hits{target};

            foreach (keys %read_hits) {
                $$stats{'reads_on_' . $_}++;
            }
        }
        $pars->close();
        $$stats{'bait_gc_hist'} = [(0) x 101 ];
        $$stats{'target_gc_hist'} = [(0) x 101];

        # add bait/target gc histogram to qc dump hash
        for my $b_or_t ('bait', 'target'){
            foreach my $chr (keys %{$intervals->{$b_or_t}}) {
                foreach my $a (@{$intervals->{$b_or_t}{$chr}}) {
                    $$stats{$b_or_t . '_gc_hist'}[$a->[2]]++;
                }
            }
        }
        $self->_update_bam_exome_qc_stats($intervals, \%pileup, $stats, $current_rname, $low_cvg);
    }


=head2 bam_exome_qc_stats

 Title   : bam_exome_qc_stats
 Usage   : $obj->bam_exome_qc_stats('in.bam', 'out.stats.dump');
 Function: Generate a stats dump file for bam files of exome data. These are used to make statistical plots of the
           data.
 Returns : boolean (true on success)
 Args    : input bam file, output filename.
           Optional input:
             MANDATORY: either all four of:
                          - bait_interval, target_interval, ref_fa, ref_fai
                        or this:
                          - load_intervals_dump_file
                        must be given in the options hash.   Info about baits/targets and reference
                        is calulcated prior to parsing the BAM file.  This info can be dumped
                        to a file to save calculating it more than once, by running first with
                        -dump_intervals, then do the QC with -load_intervals_dump_file.

             bait_interval   => filename of baits intervals file, of the form 1 bait per line, tab
                                separated:
                                reference_seq start end
                                Order doesn't matter, and baits can overlap.  This will be dealt with.
             target_interval => filename of targets intervals file, same format as bait_interval.
             ref_fa          => filename of reference fasta file
             ref_fai         => fai file of reference fasta file
             dump_intervals  => name of file to dump info to.  If this is used, the BAM is NOT QC'd.
             load_intervals_dump_file => name of file made when dump_intervals was used.  This will
                                         load all reference/bait/target info from the dump file, which
                                         is faster than creating it from scratch
             bam             => name of bam file to be QC'd.
                                MANDATORY, unless the dump_intervals option is used (in which case
                                the BAM file is ignored)
             low_cvg         => int (default 2).  Used to make the statistic low_cvg_targets.
                                Any target with all its bases' coverage less than this value will
                                be counted as low coverage.  e.g. default means that a target
                                must have coverage of at least 2 at one base to not count as a
                                low coverage target.
             near_length     => int (default 100)
                                Distance to count as 'near' to a bait or a target
             max_insert      => int (default 2000) Ignore insert sizes > this number.  Sometimes
                                (depending on the mapper) BAM files can have unrealisticly large
                                insert sizes.  This option removes the outliers, so calculated
                                mean insert size is more realistic. 
=cut

     method bam_exome_qc_stats (
     ClassName|Object $self: Str|File :$bam, Str|File :$dump, Str|File :$load_intervals_dump_file?,  Str|File :$bait_interval?, Str|File :$target_interval?, Str|File :$ref_fa?, Str|File :$ref_fai?, Bool :$name_check?) {
        my $in_bam = VRPipe::File->get(path => $bam);
        my $stats_out = VRPipe::File->get(path => $dump); 
        my $ofh = $stats_out->openw;
        $stats_out->disconnect;
        #add optional inputs later....
        my $intervals;
        my $ref_lengths; 
        #make these optional inputs:
        my $near_length = 100;# unless $opts{near_length};
#       my $low_cvg = 2;# unless $opts{low_cvg};
#       my $max_insert = 1000;# unless $opts{max_insert};
#     	my $first_sam_record = 1;
		my $counter = 0;
        my %stats = ('bait_bases', 0,
                     'bait_bases_mapped', 0,
                     'bait_near_bases_mapped', 0,
                     'bait_near_bases', 0,
                     'bases_mapped', 0,
                     'bases_mapped_reads', 0,
                     'clip_bases', 0,
                     'low_cvg_targets', 0,
                     'mapped_as_pair', 0,
                     'num_mismatches', 0,
                     'off_bait_bases', 0,
                     'raw_bases', 0,
                     'raw_reads', 0,
                     'reads_mapped', 0,
                     'reads_on_target', 0,
                     'reads_on_target_near', 0,
                     'reads_on_bait', 0,
                     'reads_on_bait_near', 0,
                     'reads_paired', 0,
                     'rmdup_reads_mapped', 0,
                     'rmdup_bases_mapped', 0,
                     'target_bases', 0,
                     'target_bases_mapped', 0,
                     'target_near_bases', 0,
                     'target_near_bases_mapped', 0);

        #my %pileup; # for the current reference sequence in the bam file ($current_rname), this
                    # will store pileup info of the target bases, near target bases, bait bases and
                    # near bait bases.  One hash for each.  reference position => number of reads
                    # covering that position.  'near' does *not* include on, i.e. a given position
                    # cannot be both near to a target and in a target (similarly for baits).

        # Need to get the target/bait info.  Either it was precomputed and dumped
        # to a file we need to load, or we need to calculate it all now.
        # (Calculating the GC of each interval is slow)
        if ($load_intervals_dump_file) {
         	$intervals = do $load_intervals_dump_file or $self->throw("Error loading data from load_intervals_dump_file $load_intervals_dump_file");
        }
    	# parse the interval files and get gc from reference fasta
    	elsif ($bait_interval and $target_interval and $ref_fa and $ref_fai) {
    	    $ref_lengths = $self->_fai_chromosome_lengths($ref_fai->path);
            $intervals->{target} = $self->($target_interval->path);
            $intervals->{target_near} = $self->_intervals_file2hash($target_interval->path, $near_length);
            $intervals->{bait} = $self->_intervals_file2hash($bait_interval);
            $intervals->{bait_near} = $self->_intervals_file2hash($bait_interval, $near_length);

            foreach my $bait_or_target qw(bait target) {
                my $pars = VRPipe::Parser->create('fasta', {file => $ref_fa});
                my $pr = $pars->parsed_record();
                while ($pars->next_record()) {
                    my $id = $pr->[0];
                    my $seq = $pr->[1];
                    exists $intervals->{$bait_or_target}{$$id} || next;
                    foreach my $interval (@{$intervals->{$bait_or_target}{$$id}}) {
                        my $length = $interval->[1] - $interval->[0] + 1;
                        my $gc_count = substr($$seq, $interval->[0] - 1, $length) =~ tr/cgCG//;
                        push @{$interval}, int(100 * $gc_count / $length);
                    } 
                }
            }
        }
        else {
            $self->throw("Error. must use either intervals_dump_file, or give all four of ref_fa, ref_fai, bait_interval and target_interval\n");
        }

        # sanity check that the reference sequence names in the intervals hash have
        # a non-empty intersection with the names in the header of the bam file.  Maybe
        # not an ideal check, but we can't expect either to be a subset of the other, so
        # intersection non-empty will have to do.
        #print $ofh "Checking names...\n";
        if ($name_check) {
			my @names;
        	my $input = `samtools view -H $bam`;
        	my @lines = split('\n', $input); 
        	for my $line_in ( @lines ) {
       			if ($line_in =~ /^\@SQ\t.*SN:(\w+)/) {
          			push @names, $1;
        		}
    		}
            my %bam_ref_names = map {$_ => 1} @names;
            my $names_ok = 0;
            foreach (keys %bam_ref_names) {
                if (exists $intervals->{target}{$_}) {
                    $names_ok = 1;
                    last;
                }
            }
            # Unfortunately, we have 'MT' in masked and unmasked versions of the same
            # genome.  In this case, the previous test doesn't work.
            # These have chromosomes named chr1, chr2, ...etc , or 1,2,... etc.
            # Check this for chromosome 1: it's a reasonable assumption that
            # there is at least one target in chromosome 1.
            if ($names_ok and ( ($bam_ref_names{1} and $intervals->{targets}{chr1}) or ($bam_ref_names{chr1} and $intervals->{targets}{1}) ) ) {
                $names_ok = 0;
            }
            $names_ok or $self->throw("Mismatch in reference sequence names in BAM file vs ref sequences containing targets."); 
        }
        # Everything is initialised.  It's time to parse the bam file to get the stats
        $self->bam_statistics($in_bam->path, \%stats, $intervals);
        # sort out other global stats
        $stats{bait_design_efficiency} = $stats{bait_bases} != 0 ? $stats{target_bases} / $stats{bait_bases} : 'NA';
        $stats{pct_mismatches} = 100 * $stats{num_mismatches} / $stats{bases_mapped_reads};
        my %qual_stats = VRPipe::Utils::Math->new()->histogram_stats($stats{insert_size_hist});
        $stats{median_insert_size} = $qual_stats{q2};
        $stats{mean_insert_size} = $qual_stats{mean};
        $stats{insert_size_sd} = $qual_stats{standard_deviation};
        $stats{mean_bait_coverage} = $stats{bait_bases_mapped} / $stats{bait_bases};
        $stats{mean_target_coverage} = $stats{target_bases_mapped} / $stats{target_bases};

        # calculate cumulative coverage of baits and targets
        foreach my $bait_or_target qw(bait target) {
            my %cov_stats = VRPipe::Utils::Math->new()->histogram_stats($stats{$bait_or_target . '_cvg_hist'});
            $stats{$bait_or_target . '_coverage_sd'} = $cov_stats{standard_deviation};

            my $base_count = 0;
            my $total_bases = $stats{$bait_or_target . '_bases'};

            foreach my $cov (reverse sort {$a <=> $b} keys %{$stats{$bait_or_target . '_cvg_hist'}}) {
                $base_count += $stats{$bait_or_target . '_cvg_hist'}{$cov};
                $stats{$bait_or_target . '_cumulative_coverage'}{$cov} = $base_count;
                $stats{$bait_or_target . '_cumulative_coverage_pct'}{$cov} = $base_count / $stats{$bait_or_target . '_bases'};
            }
        }

        # work out % bases coverage above certain values
        foreach my $bait_or_target qw(bait target) { 
            my @vals = (1, 2, 5, 10, 20, 50, 100);
            foreach (@vals) {
                $stats{$bait_or_target . '_bases_' . $_ . 'X'} = 0;
            }
            $counter++;
            foreach my $cov (sort {$a <=> $b} keys %{$stats{$bait_or_target . '_cumulative_coverage_pct'}}) {
                while ((0 < scalar @vals) and $cov >= $vals[0]) {
                    $stats{$bait_or_target . '_bases_' . $cov . 'X'} = $stats{$bait_or_target . '_cumulative_coverage_pct'}{$cov};
                    shift @vals;
                }
                last if (0 == scalar @vals);
            }
        }
        print $ofh Dumper(\%stats);
        close $ofh;
    } 

    method _fai_chromosome_lengths (ClassName|Object $self: Str|File $fai_file) {
        my $out = {};
        open(my $fh, '<', $fai_file) or $self->throw("Error opening file $fai_file: $!\n");
        while (my $line=<$fh>)
        {
            if ( !($line=~/^(\S+)\s+(\d+)\s+/) ) { $self->throw("Unexpected format of $fai_file:\n$line\n"); }
            $$out{$1} = $2;
        }
        close $fh;
        return $out;
    }
    
    method _intervals_file2hash (ClassName|Object $self: Str|File $intervals_file, Int :$near_length?) {
        my %intervals;
        my $n_length = $near_length ? $near_length : 0;
        open my $fh, $intervals_file or $self->throw("Cannot open $intervals_file: $!");
        while (my $line = <$fh>) { 
            chomp $line;
            my ($ref_id, $start, $end) = split /\t/, $line;
            $intervals{$ref_id} = () unless (exists $intervals{$ref_id});
            $start = max(0, $start - $n_length);
            $end += $n_length;
            push @{$intervals{$ref_id}}, [$start, $end];
        }
        close $fh;
        
        foreach my $ref_id (keys %intervals) {
            my $list = $intervals{$ref_id};
            @$list = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @$list;
            my $i = 0;
            while ($i < scalar @$list - 1) {
                if ($list->[$i][1] >= $list->[$i+1][0] - 1) {
                    $list->[$i + 1] = [$list->[$i][0], max($list->[$i][1], $list->[$i + 1][1])];
                    splice @$list, $i, 1;
                }
                else {
                    $i++;
                }
            }
        }
        return \%intervals;
    }
    
    
    method _update_bam_exome_qc_stats (ClassName|Object $self: Ref|HashRef $intervals, Ref|HashRef $pileup, Ref|HashRef $stats, Str $ref_id, Str $low_cvg) {
        $stats->{per_ref_sequence}{$ref_id} = {};
        my $ref_id_hash = $stats->{per_ref_sequence}{$ref_id};

        # update stats specific to the given reference sequence
        $ref_id_hash->{bait_bases} = scalar keys %{$pileup->{bait}};
        $ref_id_hash->{bait_near_bases} = scalar keys %{$pileup->{bait_near}};
        $ref_id_hash->{target_bases} = scalar keys %{$pileup->{target}};
        $ref_id_hash->{target_near_bases} = scalar keys %{$pileup->{target_near}};
        $ref_id_hash->{targets} = scalar @{$intervals->{target}->{$ref_id}};
        $ref_id_hash->{baits} = scalar @{$intervals->{bait}->{$ref_id}};

        $ref_id_hash->{bait_bases_mapped} = sum 0, values %{$pileup->{bait}};
        $ref_id_hash->{bait_near_bases_mapped} = sum 0, values %{$pileup->{bait_near}};
        $ref_id_hash->{target_bases_mapped} = sum 0, values %{$pileup->{target}};
        $ref_id_hash->{target_near_bases_mapped} = sum 0, values %{$pileup->{target_near}};

        $ref_id_hash->{target_coverage_mean} = $ref_id_hash->{target_bases_mapped} / $ref_id_hash->{target_bases};
        $ref_id_hash->{bait_coverage_mean} = $ref_id_hash->{bait_bases_mapped} / $ref_id_hash->{bait_bases};

        # update global stats (combination of all reference sequences)
        foreach (qw(bait_bases bait_near_bases target_bases target_near_bases bait_bases_mapped bait_near_bases_mapped target_bases_mapped target_near_bases_mapped targets baits)) {
            $stats->{$_} += $ref_id_hash->{$_};
        }

        # update the per base coverage data
        foreach my $bait_or_target (qw(bait target)) {
            foreach my $cov (values %{$pileup->{$bait_or_target}}) {
                $stats->{$bait_or_target . '_cvg_hist'}{$cov}++;
            }
        }

        # for each bait and target, we want to know:
        # 1) the mean coverage and the %GC
        # 2) whether the min coverage was achieved over the whole target.
        foreach my $bait_or_target qw(bait target) {
            foreach my $interval (@{$intervals->{$bait_or_target}->{$ref_id}}) {
                my $cov_sum = int( 0.5 + sum @{$pileup->{$bait_or_target}}{($interval->[0] .. $interval->[1])} );
                my $mean_cov = int(0.5 + $cov_sum / ($interval->[1] - $interval->[0] + 1) );
                my $gc = int(0.5 + $interval->[2]);        
                $stats->{'gc_vs_' . $bait_or_target . '_cvg'}->[$gc]->{$mean_cov}++;
                my $cvg_check = 1;

                foreach my $i ($interval->[0] .. $interval->[1]) {
                    if ($pileup->{$bait_or_target}{$i} >= $low_cvg) {
                        $cvg_check = 0;
                        last;
                    }                
                } 
                $stats->{'low_cvg_' . $bait_or_target . 's'}++ if $cvg_check;
            }
        }
    }
    
#     
# =head2 bam_refnames
# 
#  Title   : bam_refnames
#  Usage   : my $ref_names = $obj->bam_refnames($bam_file);
#  Function: Get the reference sequnce names from the header of a bam file
#  Returns : reference to array of names.  names will be same order as in header
#  Args    : bam filename
# 
# =cut

# 	method bam_refnames {
#     my ($self, $bam_file) = @_;
#     
#     my @names;
#     my $st = VertRes::Wrapper::samtools->new(quiet => 1, run_method => 'open');
#     my $fh = $st->view($bam_file, undef, H => 1);
#     while (<$fh>){
#         if (/^\@SQ\t.*SN:(\w+)/) {
#             push @names, $1;
#         }
#     }
# 
#     return \@names;
# }

}

1;
