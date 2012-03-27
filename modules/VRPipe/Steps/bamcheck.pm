use VRPipe::Base;

class VRPipe::Steps::bamcheck with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { bamcheck_exe => VRPipe::StepOption->get(description => 'path to your bamcheck executable',
                                                         optional => 1,
                                                         default_value => 'bamcheck'),
                 bamcheck_options => VRPipe::StepOption->get(description => 'options to bamcheck, excluding -r and -t (which are set by reference_fasta and exome_targets_file options)',
                                                             optional => 1),
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file'),
		 exome_targets_file => VRPipe::StepOption->get(description => 'absolute path to a file describing the targets/baits used for exome pulldown (tab-delimited chr,start,end, 1-based, inclusive)',
							       optional => 1)};
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', description => 'bam files', max_files => -1) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options = $self->options;
            my $bamcheck_exe = $options->{bamcheck_exe};
            my $opts = VRPipe::Steps::bamcheck->get_bamcheck_options($options);
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                my $ifile = $bam_file->path;
                my $check_file = $self->output_file(output_key => 'bamcheck_files', basename => $ifile->basename.'.bamcheck', type => 'txt');
                my $ofile = $check_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bamcheck', 'stats_from_bamcheck', ["$bamcheck_exe $opts $ifile > $ofile", $req, {output_files => [$check_file]}]);
            }
        };
    }
    method outputs_definition {
        return { bamcheck_files => VRPipe::StepIODefinition->get(type => 'txt',
                                                                 description => 'the output of bamcheck on a bam',
                                                                 max_files => -1,
                                                                 metadata => {source_bam => 'path to the bam file this bamcheck file was created from',
                                                                              lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                              bases => 'total number of base pairs',
                                                                              reads => 'total number of reads (sequences)',
                                                                              reads_mapped => 'number of reads mapped',
                                                                              bases_mapped => 'number of bases mapped',
                                                                              bases_trimmed => 'number of bases trimmed',
                                                                              forward_reads => 'number of forward reads',
                                                                              reverse_reads => 'number of reverse reads',
                                                                              reads_paired => 'number of reads paired',
                                                                              rmdup_reads => 'number of reads after duplicate removal',
                                                                              rmdup_reads_mapped => 'number of reads mapped after duplicate removal',
                                                                              rmdup_bases => 'number of bases after duplicate removal',
                                                                              rmdup_bases_mapped => 'number of bases mapped after duplicate removal',
                                                                              avg_read_length => 'the average length of reads',
                                                                              paired => '0=single ended reads only; 1=paired end reads present',
                                                                              mean_insert_size => 'mean insert size (0 if unpaired)',
                                                                              error_rate => 'the error rate',
                                                                              mean_coverage => 'the mean coverage',
                                                                              library => 'library name',
                                                                              sample => 'sample name',
                                                                              center_name => 'center name',
                                                                              platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                              study => 'name of the study, put in the DS field of the RG header line',
                                                                              optional => ['library', 'sample', 'center_name', 'platform', 'study', 'mean_insert_size']}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Creates a bamcheck file, a file of summary stats on a bam file, and also associates some of the stats as metadata on the bam file";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method get_bamcheck_options (ClassName|Object $self: HashRef $options!) {
        my $ref = file($options->{reference_fasta});
        $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
        my $opts = '-r '.$ref;
        
        my $targets = $options->{exome_targets_file};
        if ($targets) {
            $targets = file($targets);
            $self->throw("exome_targets_file must be an absolute path") unless $targets->is_absolute;
            $opts .= ' -t '.$targets;
        }
        
        my $user_opts = $options->{bamcheck_options};
        if ($user_opts) {
            if ($user_opts =~ /-r|-t/) {
                $self->throw("neither -r nor -t should be supplied as bamcheck options");
            }
            
            $opts .= ' '.$user_opts;
        }
        
        return $opts;
    }
    
    method stats_from_bamcheck (ClassName|Object $self: Str $cmd_line) {
        my ($bam_path, $check_path) = $cmd_line =~ / (\S+) > (\S+)$/;
        $bam_path || $self->throw("bad cmd line [$cmd_line]");
        my $bam_file = VRPipe::File->get(path => $bam_path);
        my $check_file = VRPipe::File->get(path => $check_path);
        
        $bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $check_file->update_stats_from_disc(retries => 3);
        if ($check_file->s) {
            my $new_meta = {};
            
            # parse the bamcheck file
            my $parser = VRPipe::Parser->create('bamcheck', {file => $check_file});
            $check_file->disconnect;
            
            if ($cmd_line =~ /-d/) {
                $new_meta->{rmdup_reads} = $parser->sequences;
                $new_meta->{rmdup_reads_mapped} = $parser->reads_mapped;
                $new_meta->{rmdup_bases} = $parser->total_length;
                $new_meta->{rmdup_bases_mapped} = $parser->bases_mapped;
                $new_meta->{rmdup_bases_mapped_c} = $parser->bases_mapped_cigar;
                $new_meta->{rmdup_bases_trimmed} = $parser->bases_trimmed;
            }
            else {
                $new_meta->{reads} = $parser->sequences;
                $new_meta->{reads_mapped} = $parser->reads_mapped;
                $new_meta->{bases} = $parser->total_length;
                $new_meta->{bases_mapped} = $parser->bases_mapped;
                $new_meta->{bases_mapped_c} = $parser->bases_mapped_cigar;
                $new_meta->{bases_trimmed} = $parser->bases_trimmed;
                
                $new_meta->{reads_paired} = $parser->reads_paired;
                $new_meta->{paired} = $parser->is_paired;
                $new_meta->{error_rate} = $parser->error_rate;
                $new_meta->{forward_reads} = $parser->first_fragments;
                $new_meta->{reverse_reads} = $parser->last_fragments;
                $new_meta->{avg_read_length} = $parser->average_length;
                $new_meta->{mean_insert_size} = $parser->insert_size_average;
                $new_meta->{sd_insert_size} = $parser->insert_size_standard_deviation;
                
                my $reads_duplicated = $parser->reads_duplicated;
                $new_meta->{rmdup_reads} = $new_meta->{reads} - $reads_duplicated;
                $new_meta->{rmdup_reads_mapped} = $new_meta->{reads_mapped} - $reads_duplicated;
                my $bases_duplicated = $parser->bases_duplicated;
                $new_meta->{rmdup_bases} = $new_meta->{bases} - $bases_duplicated;
                $new_meta->{rmdup_bases_mapped} = $new_meta->{bases_mapped} - $bases_duplicated;
                
                $new_meta->{mean_coverage} = $parser->mean_coverage;
                foreach my $cov (1, 2, 5, 10, 20, 50, 100) {
                    $new_meta->{"bases_of_${cov}X_coverage"} = $parser->cumulative_coverage($cov);
                }
            }
            
            # and get other metadata from bam header, but don't overwrite
            # existing info
            $parser = VRPipe::Parser->create('bam', {file => $bam_file});
            my %rg_info = $parser->readgroup_info();
            my @rgs = keys %rg_info;
            my $existing_meta = $bam_file->metadata;
            if (@rgs == 1) {
                my $info = $rg_info{$rgs[0]};
                unless (defined $existing_meta->{lane}) { $new_meta->{lane} = $info->{PU} ? $info->{PU} : $rgs[0]; } else { $new_meta->{lane} = $existing_meta->{lane} }
                unless (defined $existing_meta->{library}) { $new_meta->{library} = $info->{LB} if $info->{LB}; } else { $new_meta->{library} = $existing_meta->{library} }
                unless (defined $existing_meta->{sample}) { $new_meta->{sample} = $info->{SM} if $info->{SM}; } else { $new_meta->{sample} = $existing_meta->{sample} }
                unless (defined $existing_meta->{center_name}) { $new_meta->{center_name} = $info->{CN} if $info->{CN}; } else { $new_meta->{center_name} = $existing_meta->{center_name} }
                unless (defined $existing_meta->{platform}) { $new_meta->{platform} = $info->{PL} if $info->{PL}; } else { $new_meta->{platform} = $existing_meta->{platform} }
                unless (defined $existing_meta->{study}) { $new_meta->{study} = $info->{DS} if $info->{DS}; } else { $new_meta->{study} = $existing_meta->{study} }
            }
            if (@rgs != 1 || "$new_meta->{lane}" eq "1") {
                # call the name something we can be most sure is maximally
                # unique
                my $rg = $bam_path;
                $rg =~ s/\.bam$//;
                $rg =~ s/\W/_/g;
                $new_meta->{lane} = $rg;
            }
            
            $bam_file->add_metadata($new_meta);
            my %check_meta = (%{$bam_file->metadata}, source_bam => $bam_file->path->stringify);
            $check_file->add_metadata(\%check_meta);
        }
        else {
           $self->throw("$check_path failed to be made");
        }
    }
}

1;
