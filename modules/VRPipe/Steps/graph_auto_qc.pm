
=head1 NAME

VRPipe::Steps::graph_auto_qc - a step

=head1 DESCRIPTION

*** more documentation to come

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

class VRPipe::Steps::graph_auto_qc with VRPipe::StepRole {
    use VRPipe::Schema;
    
    method options_definition {
        return {
            auto_qc_gtype_regex => VRPipe::StepOption->create(
                description   => 'A regular expression to choose acceptable genotype status values',
                optional      => 1,
                default_value => '^confirmed'
            ),
            auto_qc_mapped_base_percentage => VRPipe::StepOption->create(
                description   => 'Minimum percentage of mapped bases',
                optional      => 1,
                default_value => 90
            ),
            auto_qc_duplicate_read_percentage => VRPipe::StepOption->create(
                description   => 'Maximum percentage of duplicate reads',
                optional      => 1,
                default_value => 8
            ),
            auto_qc_mapped_reads_properly_paired_percentage => VRPipe::StepOption->create(
                description   => 'Minimum percentage of the reads that are mapped which are also properly paired',
                optional      => 1,
                default_value => 80
            ),
            auto_qc_error_rate => VRPipe::StepOption->create(
                description   => 'Maximum allowed error rate',
                optional      => 1,
                default_value => 0.02
            ),
            auto_qc_overlapping_base_duplicate_percent => VRPipe::StepOption->create(
                description   => 'Maximum percent of bases duplicated due to overlapping reads of a pair',
                optional      => 1,
                default_value => 4
            ),
            auto_qc_max_ins_to_del_ratio => VRPipe::StepOption->create(
                description   => 'Maximum insert to deletion ratio',
                optional      => 1,
                default_value => '0.82'
            ),
            auto_qc_min_ins_to_del_ratio => VRPipe::StepOption->create(
                description   => 'Minimum insert to deletion ratio',
                optional      => 1,
                default_value => '0.68'
            ),
            auto_qc_max_ic_above_median => VRPipe::StepOption->create(
                description   => 'Maximimum indels per cycle, factor above median',
                optional      => 1,
                default_value => '8'
            ),
            auto_qc_genotype_min_concordance => VRPipe::StepOption->create(
                description   => 'Minimum concordance to be used when applying genotype check status',
                optional      => 1,
                default_value => 0.94
            ),
            auto_qc_genotype_min_sites => VRPipe::StepOption->create(
                description   => 'Minimum number of sites for the data to be used in the genotype analysis',
                optional      => 1,
                default_value => 0
            ),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                description => 'bam files',
                max_files   => -1,
                metadata    => { lane => 'lane name (a unique identifer for this sequencing run, aka read group)' }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            my $opts = $self->options;
            my $req  = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $file (@{ $self->inputs->{bam_files} }) {
                my $path = $file->path;
                my $lane = $file->metadata->{lane};
                
                unless ($lane) {
                    # modern steps don't store stuff in ->metadata(); check the
                    # graph db instead
                    my $schema = VRPipe::Schema->create("VRPipe");
                    my $file_graph_node = $schema->get('File', { path => $file->path->stringify });
                    if ($file_graph_node) {
                        my $lane_node = $file_graph_node->closest('VRTrack', 'Lane', direction => 'incoming');
                        if ($lane_node) {
                            $lane = $lane_node->unique;
                        }
                    }
                }
                unless ($lane) {
                    $self->throw("file " . $file->path . " lacks lane metadata");
                }
                
                my $cmd = "use VRPipe::Steps::graph_auto_qc; VRPipe::Steps::graph_auto_qc->auto_qc(bam => q[$path], lane => q[$lane], auto_qc_gtype_regex => q[$opts->{auto_qc_gtype_regex}], auto_qc_mapped_base_percentage => $opts->{auto_qc_mapped_base_percentage}, auto_qc_duplicate_read_percentage => $opts->{auto_qc_duplicate_read_percentage}, auto_qc_mapped_reads_properly_paired_percentage => $opts->{auto_qc_mapped_reads_properly_paired_percentage}, auto_qc_error_rate => $opts->{auto_qc_error_rate}, auto_qc_max_ins_to_del_ratio => $opts->{auto_qc_max_ins_to_del_ratio}, auto_qc_min_ins_to_del_ratio => $opts->{auto_qc_min_ins_to_del_ratio}, auto_qc_overlapping_base_duplicate_percent => $opts->{auto_qc_overlapping_base_duplicate_percent}, auto_qc_max_ic_above_median => $opts->{auto_qc_max_ic_above_median}, auto_qc_genotype_min_concordance => $opts->{auto_qc_genotype_min_concordance}, auto_qc_genotype_min_sites => $opts->{auto_qc_genotype_min_sites} );";
                $self->dispatch_vrpipecode($cmd, $req);
            }
        };
    }
    
    method outputs_definition {
        return {};
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Considering the summary stats for a lane in the graph db and the metadata stored on the bam file automatically decide if the lane passes the quality check and update the graph database.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method auto_qc (ClassName|Object $self: Str|File :$bam!, Str :$lane!, Str :$auto_qc_gtype_regex?, Num :$auto_qc_mapped_base_percentage?, Num :$auto_qc_duplicate_read_percentage?, Num :$auto_qc_mapped_reads_properly_paired_percentage?, Num :$auto_qc_error_rate?, Num :$auto_qc_max_ins_to_del_ratio?, Num :$auto_qc_min_ins_to_del_ratio?, Num :$auto_qc_overlapping_base_duplicate_percent?,  Num :$auto_qc_max_ic_above_median?,  Num :$auto_qc_genotype_min_concordance?,  Num :$auto_qc_genotype_min_sites? ) {
        my $bam_file = VRPipe::File->get(path => $bam);
        my $meta = $bam_file->metadata;
        
        # get the qc node from graph db
        my $schema = VRPipe::Schema->create('VRTrack');
        my $file_graph_node = $schema->get_file($bam_file->protocolless_path, $bam_file->protocol);
        unless ($file_graph_node) {
            $self->throw($bam_file->path . " was not in the graph database");
        }
        my $qc_nodes = $schema->file_qc_nodes(node => $bam_file);
        my $bam_stats = $qc_nodes->{bam_stats}->properties;
        
        my @qc_status = ();
        my ($test, $status, $reason);
        
        # check to see if the bam file contains any reads as this crashes other parts of auto_qc
        my $bam_has_seq = 1;
        if ($bam_stats->{'sequences'} == 0 && $bam_stats->{'total length'} == 0) {
            $bam_has_seq = 0;
            push @qc_status, { test => 'Empty bam file check', status => 0, reason => 'The bam file provided for this lane contains no sequences.' };
        }
        
        # we'll always fail if the npg status is failed
        my $vrtrack_metadata = $schema->vrtrack_metadata(node => $file_graph_node);
        my $npg_status = $vrtrack_metadata->{manual_qc} || $meta->{manual_qc};
        unless ($npg_status) {
            push @qc_status, { test => 'NPG QC status check', status => 0, reason => 'The lane failed the NPG QC check, so we auto-fail as well since this data will not be auto-submitted to EGA/ENA.' };
        }
        
        # genotype check results
        if (defined $auto_qc_gtype_regex) {
            # gtype_analysis metadata key must exist if the
            # bam_genotype_checking pipeline was run; otherwise,
            # if it doesn't exist, use the genotype data from graph db.
            my $gstatus;
            my $gtype_analysis = $meta->{gtype_analysis};
            if ($gtype_analysis) {
                ($gstatus) = $gtype_analysis =~ /status=(\S+) expected=(\S+) found=(\S+) (?:ratio|concordance)=(\S+)/;
            }
            elsif ($qc_nodes->{Genotype}) {
                my $genotype    = $qc_nodes->{Genotype}->properties;
                my $expected    = $genotype->{'expected_sample_name'};
                my $found       = $genotype->{'matched_sample_name'};
                my $match_count = $genotype->{'match_count'};
                my $snp_count   = $genotype->{'common_snp_count'};
                my $concordance = $match_count / $snp_count;
                if ($expected eq $found) {
                    if ($concordance > $auto_qc_genotype_min_concordance && $snp_count > $auto_qc_genotype_min_sites) {
                        $gstatus = "confirmed";
                    }
                    else {
                        $gstatus = "unconfirmed";
                    }
                }
                else {
                    $gstatus = "unknown";
                }
                my $gt_status = "status=$gstatus expected=$expected found=$found concordance=$concordance";
                $bam_file->add_metadata({ gtype_analysis => $gt_status });
            }
            if ($gstatus) {
                $status = 1;
                $reason = qq[The status is '$gstatus'.];
                if ($gstatus !~ /$auto_qc_gtype_regex/) {
                    $status = 0;
                    $reason = "The status ($gstatus) does not match the regex ($auto_qc_gtype_regex).";
                }
                push @qc_status, { test => 'Genotype check', status => $status, reason => $reason };
            }
        }
        
        # mapped bases
        if (defined $auto_qc_mapped_base_percentage && $bam_has_seq) {
            my $min = $auto_qc_mapped_base_percentage;
            $status = 1;
            my $clip_bases     = $bam_stats->{'total length'} - $bam_stats->{'bases trimmed'};
            my $bases_mapped_c = $bam_stats->{'bases mapped (cigar)'};
            my $value          = 100 * $bases_mapped_c / $clip_bases;
            $reason = sprintf "At least %.1f%% bases mapped after clipping (%.2f%%).", $min, $value;
            if ($value < $min) {
                $status = 0;
                $reason = sprintf "Less than %.1f%% bases mapped after clipping (%.2f%%).", $min, $value;
            }
            push @qc_status, { test => 'Mapped bases', status => $status, reason => $reason };
        }
        
        # duplicate reads
        if (defined $auto_qc_duplicate_read_percentage && $bam_has_seq) {
            my $max = $auto_qc_duplicate_read_percentage;
            if ($bam_stats->{'reads mapped'} && $bam_stats->{'reads mapped after rmdup'}) {
                $status = 1;
                my $mapped_reads = $bam_stats->{'reads mapped'};
                my $dup_reads    = $mapped_reads - $bam_stats->{'reads mapped after rmdup'};
                my $value        = 100 * $dup_reads / $mapped_reads;
                $reason = sprintf "%.1f%% or fewer reads were duplicates (%.2f%%).", $max, $value;
                if ($value > $max) {
                    $status = 0;
                    $reason = sprintf "More than %.1f%% reads were duplicates (%.2f%%).", $max, $value;
                }
                push @qc_status, { test => 'Duplicate reads', status => $status, reason => $reason };
            }
        }
        
        # properly paired mapped reads
        if (defined $auto_qc_mapped_reads_properly_paired_percentage && $bam_has_seq) {
            my $min = $auto_qc_mapped_reads_properly_paired_percentage;
            if ($bam_stats->{'reads mapped'} && $bam_stats->{'reads mapped after rmdup'}) {
                $status = 1;
                my $mapped_reads = $bam_stats->{'reads mapped'};
                my $paired_reads = $bam_stats->{'reads mapped after rmdup'};
                my $value        = 100 * $paired_reads / $mapped_reads;
                $reason = sprintf "%.1f%% or more reads that were mapped are in a proper pair (%.2f%%).", $min, $value;
                if ($value < $min) {
                    $status = 0;
                    $reason = sprintf "Less than %.1f%% of reads that were mapped are in a proper pair (%.2f%%).", $min, $value;
                }
                push @qc_status, { test => 'Reads mapped in a proper pair', status => $status, reason => $reason };
            }
        }
        
        # error rate
        if (defined $auto_qc_error_rate) {
            my $max = $auto_qc_error_rate;
            $status = 1;
            my $error_rate = sprintf("%0.6f", $bam_stats->{'error rate'});
            $reason = "The error rate is lower than $max ($error_rate).";
            if (defined $error_rate && $error_rate > $max) {
                $status = 0;
                $reason = "The error rate is higher than $max ($error_rate).";
            }
            push @qc_status, { test => 'Error rate', status => $status, reason => $reason };
        }
        
        # number of insertions vs deletions
        my $inum = $bam_stats->{'number of insertions'};
        my $dnum = $bam_stats->{'number of deletions'};
        if (defined $inum && defined $dnum) {
            if (defined $auto_qc_max_ins_to_del_ratio) {
                my $max = $auto_qc_max_ins_to_del_ratio;
                $status = 1;
                $reason = "The Ins/Del ratio is smaller than $max ($inum/$dnum).";
                if (!$dnum or $inum / $dnum > $max) {
                    $status = 0;
                    $reason = "The Ins/Del ratio is bigger than $max ($inum/$dnum).";
                }
                push @qc_status, { test => 'InDel ratio', status => $status, reason => $reason };
            }
            if (defined $auto_qc_min_ins_to_del_ratio) {
                my $min = $auto_qc_min_ins_to_del_ratio;
                $status = 1;
                $reason = "The Ins/Del ratio is greater than $min ($inum/$dnum).";
                if ($dnum && (!$inum or $inum / $dnum < $min)) {
                    $status = 0;
                    $reason = "The Ins/Del ratio is smaller than $min ($inum/$dnum).";
                }
                push @qc_status, { test => 'InDel ratio minimum', status => $status, reason => $reason };
            }
        }
        
        # insert size
        if (($meta->{is_paired_read} || $bam_stats->{'lane_is_paired_read'}) && $bam_has_seq) {
            $test = 'Insert size';
            
            if ($bam_stats->{'reads paired'} == 0) {
                push @qc_status, { test => $test, status => 0, reason => 'Zero paired reads, yet flagged as paired' };
            }
            elsif ($bam_stats->{'insert size average'} == 0) {
                push @qc_status, { test => $test, status => 0, reason => 'The insert size not available, yet flagged as paired' };
            }
            else {
                # 80-25 rule, check how wide the insert size distribution is and
                # whether 80% of the data lies within 25% from the max peak
                # (e.g. [mpeak*(1-0.25),mpeak*(1+0.25)])
                # only libraries can be failed based on wrong insert size. The
                # lanes are always passed as long as the insert size is
                # consistent with other lanes from the same library.
                $status = 1;
                my $amount = $bam_stats->{'inserts within 25% max peak'};
                my $range  = $bam_stats->{'peak window containing 80% of inserts'};
                if (defined $amount && defined $range) {
                    $reason = sprintf "There are %.1f%% or more inserts within %.1f%% of max peak (%.2f%%).", 80, 25, $amount;
                    if ($amount < 80) {
                        $status = 0;
                        $reason = sprintf "Fail library, less than %.1f%% of the inserts are within %.1f%% of max peak (%.2f%%).", 80, 25, $amount;
                    }
                    push @qc_status, { test => $test, status => $status, reason => $reason };
                    
                    $status = 1;
                    $reason = sprintf "%.1f%% of inserts are contained within %.1f%% of the max peak (%.2f%%).", 80, 25, $range;
                    if ($range > 25) {
                        $status = 0;
                        $reason = sprintf "Fail library, %.1f%% of inserts are not within %.1f%% of the max peak (%.2f%%).", 80, 25, $range;
                    }
                    push @qc_status, { test => 'Insert size (rev)', status => $status, reason => $reason };
                }
                else {
                    push @qc_status, { test => $test, status => 0, reason => 'The insert size not available, yet flagged as paired' };
                }
            }
        }
        
        # overlapping base duplicate percent
        # calculate the proportion of mapped bases duplicated e.g. if a fragment
        # is 160bp - then 40bp out of 200bp sequenced (or 20% of bases sequenced
        # in the fragment are duplicate sequence)
        #
        #------------->
        #          <------------
        #        160bp
        #|---------------------|
        #          |--|
        #          40bp
        if (defined $auto_qc_overlapping_base_duplicate_percent) {
            if ($bam_stats->{'dupl mapped bases'} && $bam_stats->{'total mapped bases'}) {
                my ($dup_mapped_bases, $tot_mapped_bases) = ($bam_stats->{'dupl mapped bases'}, $bam_stats->{'total mapped bases'});
                my $percent = sprintf("%0.1f", ($dup_mapped_bases * 100) / $tot_mapped_bases);
                my $max = $auto_qc_overlapping_base_duplicate_percent;
                
                $reason = "The percent of bases duplicated due to reads of a pair overlapping ($percent) is smaller than or equal to $max.";
                my $status = 1;
                if ($percent > $max) {
                    $reason = "The percent of bases duplicated due to reads of a pair overlapping ($percent) is greater than $max.";
                    $status = 0;
                }
                push @qc_status, { test => 'Overlap duplicate base percent', status => $status, reason => $reason };
            }
        }
        
        # Maximimum indels per cycle, check if the highest indels per cycle is bigger than Nx the median, where N is the parameter
        if (defined $auto_qc_max_ic_above_median) {
            $status = 1;
            $reason = "All indels per cycle less then ${auto_qc_max_ic_above_median}X of the median";
            
            # Get median and max of indel fwd/rev cycle counts
            if ($bam_stats->{'median of indel cycle counts'} && $bam_stats->{'max of indel cycle counts'}) {
                my @med = split /,/, $bam_stats->{'median of indel cycle counts'};
                my @max = split /,/, $bam_stats->{'max of indel cycle counts'};
                for (my $i = 0; $i < 4; $i++) {
                    if ($max[$i] > $auto_qc_max_ic_above_median * $med[$i]) {
                        $status = 0;
                        $reason = "Some indels per cycle exceed ${auto_qc_max_ic_above_median}X of the median";
                    }
                }
                push @qc_status, { test => 'InDels per Cycle', status => $status, reason => $reason };
            }
        }
        
        # Get overall autoqc result
        $status = 1;
        for my $stat (@qc_status) {
            if (!$stat->{status}) {
                $status = 0;
            }
        }
        $bam_file->add_metadata({ auto_qc_status => $status ? 'passed' : 'failed' }, replace_data => 1);
        
        # Add the qc status to graph db
        my $date = time();
        my %qc_stats = map { $_->{test} => ["$_->{status}", $_->{reason}] } @qc_status;
        $schema->add(
            'Auto_QC',
            {
                date => $date,
                pass => $status,
                %qc_stats
            },
            incoming => { type => 'auto_qc_status', node => $file_graph_node }
        );
    
    }

}

1;
