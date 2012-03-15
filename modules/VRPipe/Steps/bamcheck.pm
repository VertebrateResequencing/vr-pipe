use VRPipe::Base;

class VRPipe::Steps::bamcheck with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { bamcheck_exe => VRPipe::StepOption->get(description => 'path to your bamcheck executable',
                                                         optional => 1,
                                                         default_value => 'bamcheck'),
                 bamcheck_options => VRPipe::StepOption->get(description => 'options to bamcheck, excluding the value to -r which will come from reference_fasta option',
                                                             optional => 1),
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file (for GC-depth calculation, only used if bamcheck_options contains -r)',
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
                                                                              forward_reads => 'number of forward reads',
                                                                              reverse_reads => 'number of reverse reads',
                                                                              avg_read_length => 'the average length of reads',
                                                                              paired => '0=single ended reads only; 1=paired end reads present',
                                                                              mean_insert_size => 'mean insert size (0 if unpaired)',
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
        my $opts = $options->{bamcheck_options};
        if ($opts) {
            my $ref = $options->{reference_fasta};
            if ($opts =~ /-r/ && $ref) {
                unless ($opts =~ /-r $ref/) {
                    $opts =~ s/-r/-r $ref/;
                }
            }
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
            open(my $ifh, $check_path) || die "could not open $check_path\n";
            
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
            }
            
            if ($cmd_line =~ /-r/) {
                # gc-stats?
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
