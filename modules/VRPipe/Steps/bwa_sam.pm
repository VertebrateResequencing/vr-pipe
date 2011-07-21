use VRPipe::Base;

class VRPipe::Steps::bwa_sam with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file the sai files were aligned with'),
                 bwa_samse_cmd => VRPipe::StepOption->get(description => 'the near-complete bwa samse command line, including desired options, but excluding the input sai, fastq, reference, -r and -f',
                                                          optional => 1,
                                                          default_value => 'bwa samse'),
                 bwa_sampe_cmd => VRPipe::StepOption->get(description => 'the near-complete bwa sampe command line, including desired options, but excluding the input sai, fastq, reference, -r and -f; defaults to bwa sampe with -a set as appropriate for the fastq insert_size',
                                                          optional => 1) };
    }
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq',
                                                              max_files => -1,
                                                              description => 'fastq files',
                                                              metadata => {lane => 'lane name (a unique identifer for this sequencing run)',
                                                                           library => 'library name',
                                                                           sample => 'sample name',
                                                                           center_name => 'center name',
                                                                           platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                           study => 'name of the study, put in the DS field of the RG header line',
                                                                           insert_size => 'expected (mean) insert size if paired',
                                                                           reads => 'total number of reads (sequences)',
                                                                           paired => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                                                                           mate => 'if paired, the path to the fastq that is our mate',
                                                                           chunk => 'if the fastq file was produced by fastq_split Step, the chunk number',
                                                                           optional => ['mate', 'chunk', 'library', 'insert_size', 'sample', 'center_name', 'platform', 'study']}),
                 sai_files => VRPipe::StepIODefinition->get(type => 'bin',
                                                            max_files => -1,
                                                            description => 'sai files, as produced by bwa aln',
                                                            metadata => {source_fastq => 'the fastq file that was input to bwa aln to generate this sai file',
                                                                         reference_fasta => 'the reference fasta that reads were aligned to'})};
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $ref = Path::Class::File->new($self->options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my %cmds;
            my $exe;
            foreach my $ended ('se', 'pe') {
                my $cmd = $self->options->{"bwa_sam${ended}_cmd"} || 'bwa sampe';
                my @c = split(" ", $cmd);
                unless ($c[1] eq "sam$ended") {
                    $self->throw("bad bwa_sam${ended}_cmd '$cmd'");
                }
                
                if ($cmd =~ /$ref|-r|-f/) {
                    $self->throw("bwa_sam${ended}_cmd should not include the ref or -r or -f option");
                }
                
                $exe = $c[0];
                $cmds{$ended} = $cmd;
            }
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            my @fq_files = @{$self->inputs->{fastq_files}};
            my @sai_files = @{$self->inputs->{sai_files}};
            $self->throw("Mismatching number of fastq files and sai files") unless @fq_files == @sai_files;
            
            my %fqs_by_lane;
            my %fqs_by_path;
            foreach my $fq (@fq_files) {
                my $fq_meta = $fq->metadata;
                my $paired = $fq_meta->{paired};
                my $path = $fq->path->stringify;
                my $lane = $fq_meta->{lane};
                my $chunk = $fq_meta->{chunk} || 0;
                $fqs_by_path{$path} = [$lane, $chunk, $paired, $fq_meta];
                $fqs_by_lane{$lane}->{$chunk}->{$paired} = [$path];
            }
            
            my %sais_by_lane;
            foreach my $sai (@sai_files) {
                my $sai_meta = $sai->metadata;
                my $path = $sai->path->stringify;
                my $source_fastq = $sai_meta->{source_fastq};
                my $ref = $fqs_by_path{$source_fastq} || $self->throw("got a sai file $path with source_fastq $source_fastq, but that fastq was not an input fastq to bwa_sam Step");
                push(@{$fqs_by_lane{$ref->[0]}->{$ref->[1]}->{$ref->[2]}}, $path);
            }
            
            my $summary_cmd;
            while (my ($lane, $chunks) = each %fqs_by_lane) {
                while (my ($chunk, $ends) = each %$chunks) {
                    while (my ($paired, $paths) = each %$ends) {
                        next if $paired == 2;
                        
                        my @fqs = ($paths->[0]);
                        my @sais = ($paths->[1]);
                        my $fq_meta = $fqs_by_path{$fqs[0]}->[3];
                        my $reads = $fq_meta->{reads};
                        
                        my $this_cmd;
                        if ($paired == 0) {
                            $this_cmd = $cmds{se};
                            $summary_cmd ||= $this_cmd;
                        }
                        else {
                            $this_cmd = $cmds{pe};
                            my $mate_paths = $ends->{2};
                            push(@fqs, $mate_paths->[0]);
                            push(@sais, $mate_paths->[1]);
                            
                            unless ($this_cmd =~ /-a/) {
                                my $insert_size = $fq_meta->{insert_size} || 500;
                                my $max = $insert_size * 3;
                                $this_cmd .= " -a $max";
                            }
                            $summary_cmd = $this_cmd;
                            
                            $reads += $fqs_by_path{$fqs[1]}->[3]->{reads};
                        }
                        
                        # add metadata and construct RG line
                        my $rg_line = '@RG\tID:'.$lane;
                        my $sam_meta = {lane => $lane,
                                        reads => $reads,
                                        paired => $paired,
                                        mapped_fastqs => join(',', @fqs),
                                        $chunk ? (chunk => $chunk) : ()};
                        if (defined $fq_meta->{library}) {
                            my $lb = $fq_meta->{library};
                            $sam_meta->{library} = $lb;
                            $rg_line .= '\tLB:'.$lb;
                        }
                        if (defined $fq_meta->{sample}) {
                            my $sm = $fq_meta->{sample};
                            $sam_meta->{sample} = $sm;
                            $rg_line .= '\tSM:'.$sm;
                        }
                        if (defined $fq_meta->{insert_size}) {
                            my $pi = $fq_meta->{insert_size};
                            $sam_meta->{insert_size} = $pi;
                            $rg_line .= '\tPI:'.$pi;
                        }
                        if (defined $fq_meta->{center_name}) {
                            my $cn = $fq_meta->{center_name};
                            $sam_meta->{center_name} = $cn;
                            $rg_line .= '\tCN:'.$cn;
                        }
                        if (defined $fq_meta->{platform}) {
                            my $pl = $fq_meta->{platform};
                            $sam_meta->{platform} = $pl;
                            $rg_line .= '\tPL:'.$pl;
                        }
                        if (defined $fq_meta->{study}) {
                            my $ds = $fq_meta->{study};
                            $sam_meta->{study} = $ds;
                            $rg_line .= '\tDS:'.$ds;
                        }
                        
                        my $sam_file = $self->output_file(output_key => 'bwa_sam_files',
                                                          output_dir => Path::Class::File->new($fqs[0])->dir,
                                                          basename => $chunk ? "$lane.$chunk.sam" : "$lane.sam",
                                                          type => 'txt',
                                                          metadata => $sam_meta);
                        
                        $this_cmd .= " -r '$rg_line' -f ".$sam_file->path." $ref @sais @fqs";
                        $self->dispatch_wrapped_cmd('VRPipe::Steps::bwa_sam', 'sam_and_check', [$this_cmd, $req, {output_files => [$sam_file]}]);
                    }
                }
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => Path::Class::File->new($exe)->basename, version => VRPipe::StepCmdSummary->determine_version($exe, '^Version: (.+)$'), summary => $summary_cmd.' -r $rg_line -f $sam_file $reference_fasta $sai_file(s) $fastq_file(s)'));
        };
    }
    method outputs_definition {
        return { bwa_sam_files => VRPipe::StepIODefinition->get(type => 'txt',
                                                                max_files => -1,
                                                                description => 'mapped sam file(s)',
                                                                metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                             library => 'library name',
                                                                             sample => 'sample name',
                                                                             center_name => 'center name',
                                                                             platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                             study => 'name of the study, put in the DS field of the RG header line',
                                                                             insert_size => 'expected (mean) insert size if paired',
                                                                             reads => 'total number of reads (sequences)',
                                                                             paired => '0=unpaired reads were mapped; 1=paired reads were mapped',
                                                                             mapped_fastqs => 'comma separated list of the fastq file(s) that were mapped',
                                                                             chunk => 'if this was mapped with fastqs that were chunks of an original fastq, this tells you which chunk',
                                                                             optional => ['chunk', 'library', 'insert_size', 'sample', 'center_name', 'platform', 'study']}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Produces sam files with bwa samse/sampe for each single or pair of input fastqs&sais";
    }
    
    method sam_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($sam_path) = $cmd_line =~ /-f (\S+)/;
        $sam_path || $self->throw("cmd_line [$cmd_line] had no -f output specified");
        
        my $sam_file = VRPipe::File->get(path => $sam_path);
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]"); #*** want to exit with the exit code of bwa failing...
        
        my $expected_reads = $sam_file->metadata->{reads};
        $sam_file->update_stats_from_disc(retries => 3);
        my $lines = $sam_file->lines;
        
        if ($lines > $expected_reads) {
            return 1;
        }
        else {
            $sam_path->unlink;
            $self->throw("cmd [$cmd_line] failed because $lines lines were generated in the sam file, yet there were $expected_reads reads in the fastq file(s)");
        }
    }
}

1;
