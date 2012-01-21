use VRPipe::Base;

class VRPipe::Steps::bam_to_fastq extends VRPipe::Steps::picard {
    use VRPipe::Parser;
    
    around options_definition {
        return { %{$self->$orig},
                 samtofastq_options => VRPipe::StepOption->get(description => 'options for picard SamToFastq, excluding input and output specifiers', optional => 1, default_value => 'VALIDATION_STRINGENCY=SILENT')
                };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', 
                                                            max_files => -1, 
                                                            description => '1 or more bam files',
                                                            metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                         bases => 'total number of base pairs',
                                                                         reads => 'total number of reads (sequences)',
                                                                         forward_reads => 'number of forward reads',
                                                                         reverse_reads => 'number of reverse reads',
                                                                         paired => '0=single ended reads only; 1=paired end reads present',
                                                                         insert_size => 'average insert size (0 if unpaired)',
                                                                         library => 'library name',
                                                                         sample => 'sample name',
                                                                         center_name => 'center name',
                                                                         platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                         study => 'name of the study, put in the DS field of the RG header line',
                                                                         optional => ['library', 'sample', 'center_name', 'platform', 'study', 'insert_size']}
                                                            ),
                };
    }
    method body_sub {
        return sub {
            use VRPipe::Utils::picard;
            
            my $self = shift;
            my $options = $self->options;
            my $picard = VRPipe::Utils::picard->new(picard_path => $options->{picard_path}, java_exe => $options->{java_exe});
            my $stf_jar = Path::Class::File->new($picard->picard_path, 'SamToFastq.jar');
            
            my $opts = $options->{samtofastq_options};
            if ($opts =~ /SamToFastq|INPUT|FASTQ|OUTPUT_DIR/i) {
                $self->throw("bam_to_fastq_options should not include the SamToFastq task command, or input/output specifiers");
            }
            if ($opts =~ /OUTPUT_PER_RG/i) {
                $self->throw("OUTPUT_PER_RG is not supported here; you must use lane (single RG) bams");
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'picard', 
                                   version => $picard->determine_picard_version(),
                                   summary => 'java $jvm_args -jar SamToFastq.jar INPUT=$bam_file(s) FASTQ=$lane.1.fastq SECOND_END_FASTQ=$lane.2.fastq'.$opts));
            my $req = $self->new_requirements(memory => 1000, time => 1);
            
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $meta = $bam->metadata;
                my $paired = $meta->{paired};
                
                my $source_bam = $bam->path->stringify;
                my $fastq_meta = { source_bam => $source_bam };
                foreach my $key (qw(lane insert_size library sample center_name platform study)) {
                    if (defined $meta->{$key}) {
                        $fastq_meta->{$key} = $meta->{$key};
                    }
                }
                
                my $out_spec;
                my @fastqs;
                if ($paired) {
                    $fastq = $self->output_file(output_key => 'fastq_files',
                                                basename => "$fastq_meta->{lane}.1.fastq",
                                                type => 'fastq',
                                                metadata => {%$fastq_meta,
                                                             reads => $meta->{forward_reads},
                                                             paired => 1});
                    my $reverse = $self->output_file(output_key => 'fastq_files',
                                                     basename => "$fastq_meta->{lane}.2.fastq",
                                                     type => 'fastq',
                                                     metadata => {%$fastq_meta,
                                                                  reads => $meta->{reverse_reads},
                                                                  paired => 2});
                    @fastqs = ($fastq, $reverse);
                    
                    $fastq->add_metadata(mate => $reverse->path->stringify);
                    $reverse->add_metadata(mate => $fastq->path->stringify);
                    
                    $out_spec = 'FASTQ='.$fastq->path.' SECOND_END_FASTQ='.$reverse->path;
                }
                else {
                    $fastq = $self->output_file(output_key => 'fastq_files',
                                                basename => "$fastq_meta->{lane}.0.fastq",
                                                type => 'fastq',
                                                metadata => {%$fastq_meta,
                                                             reads => $meta->{reads},
                                                             bases => $meta->{bases},
                                                             avg_read_length => sprintf("%0.2f", $meta->{reads} / $meta->{bases}),
                                                             paired => 0});
                    @fastqs = ($fastq);
                    
                    $out_spec = 'FASTQ='.$fastq->path;
                }
                
                my $temp_dir = $options->{tmp_dir} || $fastq->dir;
                my $jvm_args = $picard->jvm_args($req->memory, $temp_dir);
                
                my $this_cmd = $picard->java_exe.qq[ $jvm_args -jar $stf_jar INPUT=$source_bam $out_spec $opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_to_fastq', 'bam_to_fastq_and_check', [$this_cmd, $req, {output_files => \@fastqs}]); 
            }
        };
    }
    method outputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fastq', 
                                                              max_files => -1, 
                                                              description => '1 or more fastq files',
                                                              metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                           bases => 'total number of base pairs',
                                                                           reads => 'total number of reads (sequences)',
                                                                           avg_read_length => 'the average length of reads',
                                                                           paired => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                                                                           mate => 'if paired, the path to the fastq that is our mate',
                                                                           source_bam => 'path of the bam file this fastq file was made from',
                                                                           insert_size => 'average insert size (0 if unpaired)',
                                                                           library => 'library name',
                                                                           sample => 'sample name',
                                                                           center_name => 'center name',
                                                                           platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                           study => 'name of the study',
                                                                           optional => ['mate', 'library', 'sample', 'center_name', 'platform', 'study']}
                                                                    ),
               };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Converts bam files to fastq files using Picard";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method bam_to_fastq_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($bam_path, $fastq_path) = $cmd_line =~ /INPUT=(\S+) FASTQ=(\S+)/;
        my ($reverse_path) = $cmd_line =~ /SECOND_END_FASTQ=(\S+)/;
        $fastq_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file = VRPipe::File->get(path => $bam_path);
        my @out_files;
        foreach my $fq_path ($fastq_path, $reverse_path) {
            next unless $fq_path;
            push(@out_files, VRPipe::File->get(path => $fq_path));
        }
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        # check the fastq files are as expected
        my ($actual_reads, $actual_bases) = (0, 0);
        my %extra_meta;
        foreach my $out_file (@out_files) {
            $out_file->update_stats_from_disc(retries => 3);
            
            my ($these_reads, $these_bases) = (0, 0);
            my $pars = VRPipe::Parser->create('fastq', {file => $out_file});
            my $pr = $pars->parsed_record;
            while ($pars->next_record()) {
                my $id = $pr->[0];
                my $seq_len = length($pr->[1]);
                my $qual_len = length($pr->[2]);
                unless ($seq_len == $qual_len) {
                    $out_file->unlink;
                    $self->throw("Made fastq file ".$out_file->path." but sequence $id had mismatching sequence and quality lengths ($seq_len vs $qual_len)");
                }
                $these_reads++;
                $these_bases += $seq_len;
            }
            
            my $fq_meta = $out_file->metadata;
            my $expected_reads = $fq_meta->{reads};
            unless ($expected_reads == $these_reads) {
                $self->throw("Made fastq file ".$out_file->path." but there were only $these_reads reads instead of $expected_reads");
            }
            $actual_reads += $these_reads;
            
            my $expected_bases = $fq_meta->{bases};
            if ($expected_bases) {
                unless ($expected_bases == $these_bases) {
                    $out_file->unlink;
                    $self->throw("Made fastq file ".$out_file->path." but there were only $these_bases bases instead of $expected_bases");
                }
            }
            else {
                $extra_meta{$out_file->id}->{bases} = $these_bases;
            }
            $actual_bases += $these_bases;
        }
        
        my $bam_meta = $in_file->metadata;
        unless ($actual_reads == $bam_meta->{reads}) {
            foreach my $out_file (@out_files) {
                $out_file->unlink;
            }
            $self->throw("The total reads in output fastqs was only $actual_reads instead of $bam_meta->{reads}");
        }
        unless ($actual_bases == $bam_meta->{bases}) {
            foreach my $out_file (@out_files) {
                $out_file->unlink;
            }
            $self->throw("The total bases in output fastqs was only $actual_bases instead of $bam_meta->{bases}");
        }
        
        # add extra metadata we didn't know before for paired fastqs
        foreach my $out_file (@out_files) {
            my $extra = $extra_meta{$out_file->id} || next;
            my $current_meta = $out_file->metadata;
            $out_file->add_metadata(bases => $extra->{bases},
                                    avg_read_length => sprintf("%0.2f", $current_meta->{reads} / $extra->{bases}));
        }
    }
}

1;
