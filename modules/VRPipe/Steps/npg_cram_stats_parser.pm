
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

class VRPipe::Steps::npg_cram_stats_parser with VRPipe::StepRole {
    use VRPipe::Schema;
    use VRPipe::FileProtocol;
    use VRPipe::Parser;
    use JSON::XS;
    use Digest::MD5;
    
    method options_definition {
        return {
            qc_file_suffixes => VRPipe::StepOption->create(description => 'The file name suffixes of the qc files you wish to parse, comma separated.',                                                     default_value => '_F0x900.stats,.genotype.json,.verify_bam_id.json'),
            sample_id_type   => VRPipe::StepOption->create(description => 'The type of sample identifier to check in the cram file header; allowed values are id, public_name, supplier_name or accession', default_value => 'accession')
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
                $self->dispatch_vrpipecode("use VRPipe::Steps::npg_cram_stats_parser; VRPipe::Steps::npg_cram_stats_parser->parse_stats($cram_file_id, qq[$suffixes], qq[$sample_id_type])", $req);
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
    
    method parse_stats (ClassName|Object $self: Int $cram_file_id, Str $suffixes, Str $sample_id_type) {
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
            }
            $stats_graph_node->close;
            
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
            # relevant data and store in our graph db
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
                    defined $data->{sample_name_match}
                    ? (
                        matched_sample_name => $data->{sample_name_match}->{matched_sample_name},
                        match_score         => $data->{sample_name_match}->{match_score},
                        common_snp_count    => $data->{sample_name_match}->{common_snp_count},
                        match_count         => $data->{sample_name_match}->{match_count},
                        mismatch_count      => $data->{sample_name_match}->{mismatch_count},
                      )
                    : (matched_sample_name => 'n/a')
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
        my $props              = $graph_file->properties(flatten_parents => 1);
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
        return $json->decode($str);
    }
}

1;
