
=head1 NAME

VRPipe::Steps::bwa_mem_fastq - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::bwa_mem_fastq with VRPipe::StepRole {
    method options_definition {
        return {
            reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file to map against'),
            bwa_mem_options => VRPipe::StepOption->create(
                description   => 'options to bwa mem, excluding the input fastq, and reference',
                optional      => 1,
                default_value => '-M'
            ),
            bwa_exe => VRPipe::StepOption->create(
                description   => 'path to your bwa executable',
                optional      => 1,
                default_value => 'bwa'
            )
        };
    }
    
    method inputs_definition {
        return {
            fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => -1,
                description => 'fastq files, which will be alnd independently',
                metadata    => {
                    lane           => 'lane name (a unique identifer for this sequencing run)',
                    library        => 'library name',
                    sample         => 'sample name',
                    center_name    => 'center name',
                    platform       => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study          => 'name of the study, put in the DS field of the RG header line',
                    insert_size    => 'expected (mean) insert size if paired',
                    analysis_group => 'project analysis group',
                    population     => 'sample population',
                    bases          => 'total number of base pairs',
                    reads          => 'total number of reads (sequences)',
                    paired         => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                    mate           => 'if paired, the path to the fastq that is our mate',
                    chunk          => 'if the fastq file was produced by fastq_split Step, the chunk number',
                    optional       => ['mate', 'chunk', 'library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study']
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            my $ref     = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $bwa_exe  = $options->{bwa_exe};
            my $bwa_opts = $options->{bwa_mem_options};
            if ($bwa_opts =~ /$ref|mem/) {
                $self->throw("bwa_mem_options should not include the reference or the mem sub command");
            }
            my $cmd = $bwa_exe . ' mem ' . $bwa_opts;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'bwa', version => VRPipe::StepCmdSummary->determine_version($bwa_exe, '^Version: (.+)$'), summary => 'bwa mem ' . $bwa_opts . ' $reference_fasta $fastq_file(s) > $sam_file'));
            
            my ($cpus) = $bwa_opts =~ m/-t\s*(\d+)/;
            my $req = $self->new_requirements(memory => 4900, time => 2, $cpus ? (cpus => $cpus) : ());
            
            my @fq_files = @{ $self->inputs->{fastq_files} };
            
            my %fqs_by_lane;
            my %fqs_by_path;
            foreach my $fq (@fq_files) {
                my $fq_meta = $fq->metadata;
                my $paired  = $fq_meta->{paired};
                my $path    = $fq->resolve->path->stringify; # since we're storing this path in output metadata, we want the real path of the fq, not that of a symlink
                my $lane    = $fq_meta->{lane};
                my $chunk   = $fq_meta->{chunk} || 0;
                $fqs_by_path{$path} = [$lane, $chunk, $paired, $fq_meta];
                $fqs_by_lane{$lane}->{$chunk}->{$paired} = $path;
            }
            
            while (my ($lane, $chunks) = each %fqs_by_lane) {
                while (my ($chunk, $ends) = each %$chunks) {
                    while (my ($paired, $path) = each %$ends) {
                        next if $paired == 2;
                        
                        my @fqs     = ($path);
                        my $fq_meta = $fqs_by_path{ $fqs[0] }->[3];
                        my $reads   = $fq_meta->{reads};
                        my $bases   = $fq_meta->{bases};
                        
                        if ($paired != 0) {
                            my $mate_path = $ends->{2};
                            push(@fqs, $mate_path);
                            $reads += $fqs_by_path{ $fqs[1] }->[3]->{reads};
                            $bases += $fqs_by_path{ $fqs[1] }->[3]->{bases};
                        }
                        
                        # add metadata and construct RG line
                        my $rg_line  = '@RG\tID:' . $lane;
                        my $sam_meta = {
                            lane          => $lane,
                            bases         => $bases,
                            reads         => $reads,
                            paired        => $paired,
                            mapped_fastqs => join(',', @fqs),
                            $chunk ? (chunk => $chunk) : ()
                        };
                        if (defined $fq_meta->{library}) {
                            my $lb = $fq_meta->{library};
                            $sam_meta->{library} = $lb;
                            $rg_line .= '\tLB:' . $lb;
                        }
                        if (defined $fq_meta->{sample}) {
                            my $sm = $fq_meta->{sample};
                            $sam_meta->{sample} = $sm;
                            $rg_line .= '\tSM:' . $sm;
                        }
                        if (defined $fq_meta->{insert_size}) {
                            my $pi = $fq_meta->{insert_size};
                            $sam_meta->{insert_size} = $pi;
                            $rg_line .= '\tPI:' . $pi;
                        }
                        if (defined $fq_meta->{center_name}) {
                            my $cn = $fq_meta->{center_name};
                            $sam_meta->{center_name} = $cn;
                            $rg_line .= '\tCN:' . $cn;
                        }
                        if (defined $fq_meta->{platform}) {
                            my $pl = $fq_meta->{platform};
                            $sam_meta->{platform} = $pl;
                            $rg_line .= '\tPL:' . $pl;
                        }
                        if (defined $fq_meta->{study}) {
                            my $ds = $fq_meta->{study};
                            $sam_meta->{study} = $ds;
                            $rg_line .= '\tDS:' . $ds;
                        }
                        if (defined $fq_meta->{analysis_group}) {
                            $sam_meta->{analysis_group} = $fq_meta->{analysis_group};
                        }
                        if (defined $fq_meta->{population}) {
                            $sam_meta->{population} = $fq_meta->{population};
                        }
                        
                        my $ended = $paired ? 'pe' : 'se';
                        my $sam_file = $self->output_file(
                            output_key => 'bwa_mem_sam_files',
                            basename   => $chunk ? "$lane.$ended.$chunk.sam" : "$lane.$ended.sam",
                            type       => 'txt',
                            metadata   => $sam_meta
                        );
                        
                        my $this_cmd = "$cmd -R '$rg_line' $ref @fqs > " . $sam_file->path;
                        $self->dispatch_wrapped_cmd('VRPipe::Steps::bwa_mem_fastq', 'mem_and_check', [$this_cmd, $req, { output_files => [$sam_file] }]);
                    }
                }
            }
        };
    }
    
    method outputs_definition {
        return {
            bwa_mem_sam_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'output sam file(s) mapped by bwa mem',
                metadata    => {
                    lane           => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    library        => 'library name',
                    sample         => 'sample name',
                    center_name    => 'center name',
                    platform       => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study          => 'name of the study, put in the DS field of the RG header line',
                    insert_size    => 'expected (mean) insert size if paired',
                    analysis_group => 'project analysis group',
                    population     => 'sample population',
                    bases          => 'total number of base pairs',
                    reads          => 'total number of reads (sequences)',
                    paired         => '0=unpaired reads were mapped; 1=paired reads were mapped',
                    mapped_fastqs  => 'comma separated list of the fastq file(s) that were mapped',
                    chunk          => 'if this was mapped with fastqs that were chunks of an original fastq, this tells you which chunk',
                    optional       => ['chunk', 'library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study']
                }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Aligns the input fastq(s) with bwa mem to the reference";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method mem_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($sam_path) = $cmd_line =~ /> (\S+)$/;
        $sam_path || $self->throw("cmd_line [$cmd_line] was not contructed as expected");
        
        my $sam_file = VRPipe::File->get(path => $sam_path);
        
        $sam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $expected_reads = $sam_file->metadata->{reads};
        $sam_file->update_stats_from_disc(retries => 3);
        
        my $actual_records = 0;
        my $ft = VRPipe::FileType->create('bam', { file => $sam_path });
        $sam_file->disconnect;
        $actual_records = $ft->num_records;
        
        if ($actual_records == $expected_reads) {
            return 1;
        }
        elsif ($actual_records > $expected_reads) {
            # allow multiple alignments
            $sam_file->add_metadata({ reads => $actual_records }, replace_data => 1);
            return 1;
        }
        else {
            $sam_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_records records were generated in the sam file, yet there were $expected_reads reads in the fastq file(s)");
        }
    }
}

1;
