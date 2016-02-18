
=head1 NAME

VRPipe::Steps::npg_cram_stats_parser - a step

=head1 DESCRIPTION

This step is Sanger-specific.

Input should be cram files from an irods datasource with no local_root_dir
option set. It will then parse the associated stats files.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::npg_cram_stats_parser  extends VRPipe::Steps::irods {
    use VRPipe::Schema;
    use VRPipe::FileProtocol;
    use VRPipe::Parser;
    use JSON::XS;
    use Digest::MD5;
    use VRPipe::Persistent::InMemory;
    
    method _build_irods_exes {
        return { iget => 'iget', ichksum => 'ichksum', imeta => 'imeta' };
    }
    
    around options_definition {
        return {
            %{ $self->$orig },
            irods_get_zone => VRPipe::StepOption->create(
                description   => 'the zone (top level directory) where your official genotyping snp file is stored in iRODs',
                optional      => 1,
                default_value => 'seq'
            ),
            qc_file_suffixes                => VRPipe::StepOption->create(description => 'The file name suffixes of the qc files you wish to parse, comma separated.',                                                                                                               default_value => '_F0xB00.stats,.genotype.json,.verify_bam_id.json'),
            sample_id_type                  => VRPipe::StepOption->create(description => 'The type of sample identifier to check in the cram file header; allowed values are id, public_name, supplier_name or accession',                                                           default_value => 'accession'),
            genotyping_snp_file             => VRPipe::StepOption->create(description => 'Absolute path to a local file containing the correctly ordered and complemented alleles that are checked by NPG\'s genotype checker (optional if an official file is available in irods)', optional      => 1),
            genotyping_snp_file_storage_dir => VRPipe::StepOption->create(description => 'absolute path to a directory where the official genotyping snp file can be downloaded to',                                                                                                 optional      => 1),
        };
    }
    
    method inputs_definition {
        return { cram_files => VRPipe::StepIODefinition->create(type => 'cram', description => 'cram files', max_files => -1) };
    }
    
    method body_sub {
        return sub {
            my $self           = shift;
            my $opts           = $self->options;
            my $sample_id_type = $opts->{sample_id_type};
            my $suffixes       = $opts->{qc_file_suffixes};
            my $user_snp_file  = $opts->{genotyping_snp_file} || '';
            my $manifest_dir   = $opts->{genotyping_snp_file_storage_dir} || '';
            my $imeta          = $opts->{imeta_exe};
            my $iget           = $opts->{iget_exe};
            my $ichksum        = $opts->{ichksum_exe};
            my $zone           = $opts->{irods_get_zone};
            my $schema         = VRPipe::Schema->create('VRTrack');
            
            unless ($sample_id_type =~ /^(?:id|public_name|supplier_name|accession)$/) {
                $self->throw("Invalid sample_id_type $sample_id_type; must be one of id|public_name|supplier_name|accession");
            }
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $stats_suffix;
            foreach my $suffix (split(/,/, $suffixes)) {
                if ($suffix =~ /stats/) {
                    $stats_suffix = $suffix;
                    last;
                }
            }
            
            foreach my $cram_file (@{ $self->inputs->{cram_files} }) {
                # (breaking step independence) we're going to claim that the
                # stats file associated with $cram_file is our output file, so
                # that the plot_bamstats step later in our pipeline has an
                # input to work on
                my $graph_file = $schema->get_file($cram_file->protocolless_path, $cram_file->protocol);
                $self->throw($cram_file->path . " was not in the graph database") unless $graph_file;
                
                if ($stats_suffix) {
                    my $stats_graph_node;
                    foreach my $qc_file ($graph_file->related(outgoing => { type => 'qc_file' })) {
                        my $qc_path = $qc_file->path;
                        if ($qc_path =~ /^\S+$stats_suffix$/) {
                            $stats_graph_node = $qc_file;
                        }
                    }
                    $self->throw("No qc file matching *$stats_suffix was found associated with " . $graph_file->path) unless $stats_graph_node;
                    my $f = file($stats_graph_node->protocolless_path);
                    $self->output_file(output_key => 'stats_files', basename => $f->basename, output_dir => $f->dir, protocol => $stats_graph_node->protocol, type => 'txt');
                }
                
                my $cram_file_id = $cram_file->id;
                $self->dispatch_vrpipecode("use VRPipe::Steps::npg_cram_stats_parser; VRPipe::Steps::npg_cram_stats_parser->parse_stats($cram_file_id, qq[$suffixes], qq[$sample_id_type], qq[$user_snp_file], qq[$manifest_dir], qq[$imeta], qq[$iget], qq[$ichksum], qq[$zone])", $req);
            }
        };
    }
    
    method outputs_definition {
        return {
            stats_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'the samtools stats file that was associated with the input cram',
                max_files   => -1,
                min_files   => 0
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Parse qc files that have been associated with the input cram file by the irods datasource, storing the information in the graph db. Also compares the cram header to stored information and records a diff of essential data. NB: HTSLIB environment variable must be set.";
    }
    
    method max_simultaneous {
        return 50;
    }
    
    method parse_stats (ClassName|Object $self: Int $cram_file_id, Str $suffixes, Str $sample_id_type, Str $user_snp_file, Str $manifest_dir, Str $imeta, Str $iget, Str $ichksum, Str $zone) {
        my $cram_file = VRPipe::File->get(id => $cram_file_id);
        my @desired_suffixes = split(/,/, $suffixes) if $suffixes;
        my $desired_suffixes = join('|', @desired_suffixes) if @desired_suffixes;
        my $desired_qc_file_regex = qr/\S+(?:$desired_suffixes)/ if $desired_suffixes;
        
        my $schema = VRPipe::Schema->create('VRTrack');
        my $graph_file = $schema->get_file($cram_file->protocolless_path, $cram_file->protocol);
        $self->throw($cram_file->path . " was not in the graph database") unless $graph_file;
        
        # we know how to parse 3 kinds of file; pick these out based on user's
        # desired suffixes
        my ($stats_graph_node, $geno_graph_node, $verify_graph_node);
        foreach my $qc_file ($graph_file->related(outgoing => { type => 'qc_file' })) {
            my $qc_path = $qc_file->path;
            next unless $qc_path =~ /^$desired_qc_file_regex/;
            
            if ($qc_path =~ /\.stats$/) {
                if ($stats_graph_node) {
                    $self->throw("$qc_path matched $desired_suffixes and seems to be a stats file, but we already matched " . $stats_graph_node->path);
                }
                $stats_graph_node = $qc_file;
            }
            elsif ($qc_path =~ /\.genotype\.json$/) {
                if ($geno_graph_node) {
                    $self->throw("$qc_path matched $desired_suffixes and seems to be a genotype file, but we already matched " . $geno_graph_node->path);
                }
                $geno_graph_node = $qc_file;
            }
            elsif ($qc_path =~ /\.verify_bam_id\.json$/) {
                if ($verify_graph_node) {
                    $self->throw("$qc_path matched $desired_suffixes and seems to be a verify_bam_id file, but we already matched " . $verify_graph_node->path);
                }
                $verify_graph_node = $qc_file;
            }
            else {
                $self->throw("$qc_path matched $desired_suffixes but doesn't seem to be a supported type of file we can parse");
            }
        }
        
        my $date = time();
        if ($stats_graph_node) {
            # parse the stats file (don't use our bamcheck parser so we just
            # store whatever new or changed field is in the stats file)
            my $fh    = $stats_graph_node->openr;
            my %stats = ();
            my $mode  = 'normal';
            my $opts;
            my $insertion_counts = 0;
            my $deletion_counts  = 0;
            my @indels_cycles;
            my @insert_size;
            my ($short_paired_reads, $normal_paired_reads, $total_paired_reads, $dup_mapped_bases, $tot_mapped_bases);
            
            while (<$fh>) {
                if (/^# The command line was:\s+stats\s+(.+)/) {
                    $opts = $1;
                    $opts =~ s/\s+-\s*$//;
                    $opts =~ s/^\s+//;
                    $opts =~ s/\s+$//;
                    if ($opts =~ /(?<!\S)(?:-d|--remove-dups)(?!\S)/) {
                        $mode = 'rmdup';
                    }
                    elsif ($opts =~ /(?<!\S)-t(?!\S)/) {
                        $mode = 'targeted';
                    }
                }
                elsif (/^SN\s+([^:]+):\s+(\S+)/) {
                    $stats{$1} = $2;
                    
                    if ($1 eq 'raw total sequences' && $mode eq 'targeted') {
                        # we don't know the true total number of sequences
                        # in the whole file, only in the targeted regions
                        $stats{'targeted total sequences'} = delete $stats{'raw total sequences'};
                    }
                }
                elsif (/^ID\s+(\S+)\s+(\S+)\s+(\S+)/) {
                    $insertion_counts += $2;
                    $deletion_counts  += $3;
                }
                elsif (/^IC\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
                    my @items = ($2, $3, $4, $5);
                    push @indels_cycles, \@items;
                }
                elsif (/^IS\s+(\S+)\s+(\S+)/) {
                    my ($is, $pairs_total) = ($1, $2);
                    my @items = ($1, $2);
                    push @insert_size, \@items;
                    next unless $pairs_total;
                    my $seqlen = $stats{'average length'};
                    if (($seqlen * 2) > $is) {
                        $short_paired_reads += $pairs_total;
                        $dup_mapped_bases += $pairs_total * (($seqlen * 2) - $is);
                    }
                    else {
                        $normal_paired_reads += $pairs_total;
                    }
                    $total_paired_reads += $pairs_total;
                    $tot_mapped_bases += $pairs_total * ($seqlen * 2);
                }
            }
            $stats_graph_node->close;
            
            $stats{'number of insertions'} = $insertion_counts;
            $stats{'number of deletions'}  = $deletion_counts;
            
            # Calculate median and max of indel fwd/rev cycle counts
            if (@indels_cycles) {
                my (@med, @max);
                for (my $i = 0; $i < 4; $i++) {
                    my @counts = map  { $$_[$i] } @indels_cycles;
                    my @sorted = sort { $a <=> $b } @counts;
                    my $n      = int(scalar @sorted / 2);
                    push(@med, $sorted[$n]);
                    push(@max, $sorted[-1]);
                }
                $stats{'median of indel cycle counts'} = join(",", @med);
                $stats{'max of indel cycle counts'}    = join(",", @max);
            }
            
            $stats{'dupl mapped bases'}  = $dup_mapped_bases;
            $stats{'total mapped bases'} = $tot_mapped_bases;
            
            # for quick QC here we pre-compute the number of inserts that
            # fall within 25% of maximum insert peak and the window around
            # the max peak that contains 80% of the inserts
            if (@insert_size) {
                my ($amount, $range) = $self->insert_size_allowed_amount_and_range(\@insert_size);
                $stats{'inserts within 25% max peak'}           = $amount;
                $stats{'peak window containing 80% of inserts'} = $range;
            }
            
            unless ($mode eq 'rmdup') {
                # when not using -d we can still add some rmdup stats
                my $reads_duplicated = $stats{'reads duplicated'};
                $stats{'reads after rmdup'}        = $stats{'raw total sequences'} - $reads_duplicated;
                $stats{'reads mapped after rmdup'} = $stats{'reads mapped'} - $reads_duplicated;
                my $bases_duplicated = $stats{'bases duplicated'};
                $stats{'bases after rmdup'} = $stats{'total length'} - $bases_duplicated;
                
                # Estimate rmdup_bases_mapped CIGAR (excluding soft clipped
                # and for exomes parts of the reads which fall outside of
                # the target regions ) Calculate fraction
                # bases_mapped/bases_duplicated and apply this to
                # bases_mapped_cigar
                my $bases_mapped_cigar = $stats{'bases mapped (cigar)'};
                my $bases_dup_fraction = $stats{'bases mapped'} ? $bases_duplicated / $stats{'bases mapped'} : 0;
                my $bases_dup_estimate = int($bases_mapped_cigar * $bases_dup_fraction);
                $stats{'bases mapped after rmdup'} = $bases_mapped_cigar - $bases_dup_estimate;
                
                # coverage stats we just don't bother having if in -d mode;
                # we use the bamcheck parser for this since it's non-trivial
                $fh = $stats_graph_node->openr;
                my $parser = VRPipe::Parser->create('bamcheck', { file => $fh });
                $stats{'mean coverage'} = $parser->mean_coverage;
                foreach my $cov (1, 2, 5, 10, 20, 50, 100) {
                    $stats{"bases of ${cov}X coverage"} = $parser->cumulative_coverage($cov);
                }
                $stats_graph_node->close;
            }
            
            # store the parsed stats in a Bam_Stats node attached to the
            # file; because they are given an automatic new uuid as the
            # primary identifier for the node, we always end up adding a
            # new bam_stats node instead of updating an existing one. So we
            # store the current date so the user can pick the most recent
            # one
            $schema->add(
                'Bam_Stats',
                {
                    mode    => $mode,
                    options => $opts,
                    date    => $date,
                    %stats
                },
                incoming => { type => 'summary_stats', node => $stats_graph_node }
            );
        }
        
        if ($geno_graph_node) {
            # parse this json file
            my $data = $self->json_node_to_data($geno_graph_node);
            
            # based on what NPG display on their qc website, pull out the
            # relevant data and store in our graph db. The critical
            # 'Sample Reception calls' data can only be calculated based on
            # some arbitrary data stored in a file in irods, which we try to
            # get now
            my $reception_call_string = 'N' x 52;
            my $call_set              = $data->{snp_call_set};
            if (defined $data->{sample_name_match}->{matched_gt}) {
                my @alleles;
                my $parse_user_snp_file = 0;
                if (0 && $call_set && $manifest_dir) { #*** disabled while the official snp file in irods does not actually exist
                    my $snp_manifest_basename = $call_set . '.snp.manifest';
                    my $snp_manifest = VRPipe::File->create(path => file($manifest_dir, $snp_manifest_basename));
                    
                    unless ($snp_manifest->s) {
                        my ($irods_path) = $self->search_by_metadata(metadata => { snp_call_set => $call_set }, imeta => $imeta, zone => $zone);
                        if ($irods_path) {
                            # we could have multiple processes in parallel all
                            # trying to get this same file at the same time; lock
                            # and block
                            my $im       = VRPipe::Persistent::InMemory->new;
                            my $lock_key = 'npg_cram_stats_parser.' . $snp_manifest->path;
                            $im->block_until_locked($lock_key);
                            $snp_manifest->reselect_values_from_db;
                            $snp_manifest->update_stats_from_disc;
                            unless ($snp_manifest->s) {
                                $self->get_file(source => $irods_path, dest => $snp_manifest->path, iget => $iget, ichksum => $ichksum);
                            }
                            $im->unlock($lock_key);
                        }
                        elsif ($user_snp_file && -s $user_snp_file) {
                            warn "Could not find the official snp file in zone $zone matching metadata snp_call_set => $call_set, will use the user snp file\n";
                            $parse_user_snp_file = 1;
                        }
                        else {
                            $self->throw("Could not find the official snp file in zone $zone matching metadata snp_call_set => $call_set, and no user snp file was supplied");
                        }
                    }
                    
                    unless ($parse_user_snp_file) {
                        # parse the snp manifest file in to our alleles
                        my $fh = $snp_manifest->openr;
                        
                        $snp_manifest->close;
                        $self->throw("Parsing of the official snp manifest file not yet implemented, though it was downloaded to " . $snp_manifest->path);
                    }
                }
                elsif ($user_snp_file && -s $user_snp_file) {
                    $parse_user_snp_file = 1;
                }
                
                if ($parse_user_snp_file) {
                    open(my $fh, '<', $user_snp_file) || die "Failed to open $user_snp_file\n";
                    chomp(@alleles = <$fh>);
                    close $fh;
                }
                
                if (@alleles) {
                    # Each digit of the "matched_gt" is a call:
                    # 0 - NN
                    # 1 - HomA
                    # 2 - HomB
                    # 3 - Het
                    my @calls = split('', $data->{sample_name_match}->{matched_gt});
                    $self->throw("There were " . scalar(@calls) . " calls, but the snp manifest has " . scalar(@alleles) . " alleles") unless @calls == @alleles;
                    
                    $reception_call_string = '';
                    for my $i (0 .. $#alleles) {
                        my $allele = $alleles[$i];
                        my ($aa, $ab) = split('', $allele);
                        my $call = $calls[$i];
                        if ($call == 0) {
                            $reception_call_string .= '--';
                        }
                        elsif ($call == 1) {
                            $reception_call_string .= "$aa$aa";
                        }
                        elsif ($call == 2) {
                            $reception_call_string .= "$ab$ab";
                        }
                        elsif ($call == 3) {
                            $reception_call_string .= $allele;
                        }
                    }
                }
                else {
                    die "Failed to discover the alleles used for genotyping\n";
                }
            }
            
            $schema->add(
                'Genotype',
                {
                    date                     => $date,
                    expected_sample_name     => $data->{expected_sample_name},
                    pass                     => $data->{pass} ? $data->{pass} : 0,
                    bam_gt_likelihood_string => $data->{bam_gt_likelihood_string},
                    bam_gt_depths_string     => $data->{bam_gt_depths_string},
                    bam_call_count           => $data->{bam_call_count},
                    bam_call_string          => $data->{bam_call_string},
                    reception_call_string    => $reception_call_string,
                    defined $data->{sample_name_match}
                    ? (
                        matched_sample_name => $data->{sample_name_match}->{matched_sample_name},
                        match_score         => $data->{sample_name_match}->{match_score},
                        common_snp_count    => $data->{sample_name_match}->{common_snp_count},
                        match_count         => $data->{sample_name_match}->{match_count},
                        mismatch_count      => $data->{sample_name_match}->{mismatch_count},
                        matched_gt          => $data->{sample_name_match}->{matched_gt}
                      )
                    : (matched_sample_name => 'n/a'),
                    snp_call_set => $call_set
                },
                incoming => { type => 'genotype_data', node => $geno_graph_node }
            );
        }
        
        if ($verify_graph_node) {
            # parse this json file
            my $data = $self->json_node_to_data($verify_graph_node);
            
            # based on what NPG display on their qc website, pull out the
            # relevant data and store in our graph db
            $schema->add(
                'Verify_Bam_ID',
                {
                    date           => $date,
                    freemix        => $data->{freemix},
                    pass           => $data->{pass} ? $data->{pass} : 0,
                    number_of_snps => $data->{number_of_snps},
                    avg_depth      => $data->{avg_depth},
                    freeLK0        => $data->{freeLK0},
                    freeLK1        => $data->{freeLK1}
                },
                incoming => { type => 'verify_bam_id_data', node => $verify_graph_node }
            );
        }
        
        # compare the cram header to what we know about the sequencing
        my $header_lines       = $cram_file->header_lines;                                                                                                       # automagically works with irods files if HTSLIB has been compiled with irods support
        my $props              = $schema->node_and_hierarchy_properties($graph_file);
        my %rg_key_to_prop_key = (LB => ['vrtrack_library_id', 'vrtrack_library_name'], SM => 'vrtrack_sample_accession', DS => "vrtrack_study_$sample_id_type");
        my (%diffs, $ref_md5s);
        foreach (@$header_lines) {
            if (/^\@RG/) {
                while (my ($rg_key, $prop_keys) = each %rg_key_to_prop_key) {
                    my $ok = 0;
                    foreach my $prop_key (ref($prop_keys) ? @$prop_keys : ($prop_keys)) {
                        if (/\t$rg_key:([^\t]+)/ && defined $props->{$prop_key}) {
                            my $val = $1;
                            
                            if ($rg_key eq 'DS') {
                                ($val) = $val =~ /^Study ([^:]+)/;
                                $val || next;
                            }
                            
                            if ("$val" ne "$props->{$prop_key}") {
                                # like a test, we store [actual, expected]
                                $diffs{$rg_key} = [$val, $props->{$prop_key}] unless exists $diffs{$rg_key};
                            }
                            else {
                                $ok = 1;
                            }
                        }
                    }
                    
                    delete $diffs{$rg_key} if $ok;
                }
            }
            elsif (/^\@SQ.*\tM5:(\S+)/) {
                # even though we know the reference it was supposed to be
                # aligned to, the absolute path of that might differ to the abs
                # path in the SQ line, so instead we'll store an md5 of all the
                # md5s of the sequences here, and later some auto qc step could
                # check all the sequence-md5s for crams in a study is the same
                # (as each other, or the expected fasta)
                $ref_md5s .= $1;
            }
        }
        if ($ref_md5s) {
            my $dmd5 = Digest::MD5->new();
            $dmd5->add($ref_md5s);
            $ref_md5s = $dmd5->hexdigest;
        }
        $schema->add(
            'Header_Mistakes',
            {
                num_mistakes => scalar keys %diffs,
                $ref_md5s ? (md5_of_ref_seq_md5s => $ref_md5s) : (),
                %diffs
            },
            incoming => { type => 'header_mistakes', node => $graph_file }
        );
    }
    
    method json_node_to_data (ClassName|Object $self: Object $node) {
        my $json = JSON::XS->new;
        my $fh   = $node->openr;
        my $str  = do { local $/; <$fh> };
        $node->close;
        unless ($str) {
            die "Could not read any content from ", $node->path, "\n";
        }
        return $json->decode($str);
    }
    
    method insert_size_allowed_amount_and_range (ClassName|Object $self: ArrayRef $vals, Num $maxpeak_range = 25, Num $data_amount = 80) {
        # determine the max peak
        my $count       = 0;
        my $imaxpeak    = 0;
        my $ndata       = scalar @$vals;
        my $total_count = 0;
        my $max         = 0;
        for (my $i = 1; $i < $ndata; $i++) {   # skip IS=0
            my $xval = $$vals[$i][0];
            my $yval = $$vals[$i][1];
            
            $total_count += $yval;
            if ($max < $yval) {
                $imaxpeak = $i;
                $max      = $yval;
            }
        }
        
        return unless $total_count;
        
        # see how many reads are within the max peak range
        $maxpeak_range *= 0.01;
        $count = 0;
        for (my $i = 1; $i < $ndata; $i++) { # skip IS=0
            my $xval = $$vals[$i][0];
            my $yval = $$vals[$i][1];
            
            if ($xval < $$vals[$imaxpeak][0] * (1 - $maxpeak_range)) { next; }
            if ($xval > $$vals[$imaxpeak][0] * (1 + $maxpeak_range)) { next; }
            $count += $yval;
        }
        my $out_amount = 100 * $count / $total_count;
        
        # how big must be the range in order to accommodate the requested amount
        # of data
        $data_amount *= 0.01;
        my $idiff = 0;
        $count = $$vals[$imaxpeak][1];
        while ($count / $total_count < $data_amount) {
            $idiff++;
            if ($idiff < $imaxpeak)          { $count += $$vals[$imaxpeak - $idiff][1]; } # skip IS=0
            if ($idiff + $imaxpeak < $ndata) { $count += $$vals[$imaxpeak + $idiff][1]; }
            
            # this should never happen, unless $data_range is bigger than 100%
            if ($idiff > $imaxpeak && $idiff + $imaxpeak >= $ndata) { last; }
        }
        my $out_range  = $idiff <= $imaxpeak         ? $$vals[$imaxpeak][0] - $$vals[$imaxpeak - $idiff][0] : $$vals[$imaxpeak][0];
        my $out_range2 = $idiff + $imaxpeak < $ndata ? $$vals[$imaxpeak + $idiff][0] - $$vals[$imaxpeak][0] : $$vals[-1][0] - $$vals[$imaxpeak][0];
        if ($out_range2 > $out_range) { $out_range = $out_range2; }
        $out_range = 100 * $out_range / $$vals[$imaxpeak][0];
        return ($out_amount, $out_range);
    }
}

1;
