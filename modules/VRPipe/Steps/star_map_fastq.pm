
=head1 NAME

VRPipe::Steps::star_map_fastq - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014-15 Genome Research Limited.

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

class VRPipe::Steps::star_map_fastq with VRPipe::StepRole {
    method options_definition {
        return {
            star_map_options => VRPipe::StepOption->create(
                description   => 'options for star mapping, excluding input fastq(s) and genome reference/dir options',
                optional      => 1,
                default_value => '--outSAMunmapped Within --runThreadN 16'
            ),
            star_genomeGenerate_options => VRPipe::StepOption->create(
                description   => 'options to STAR genomeGenerate, excluding genome directory option (default specifies parameters required to estimate the memory size)',
                optional      => 1,
                default_value => '--genomeSAindexNbases 14 --runThreadN 16 --limitIObufferSize 150000000'
            ),
            star_exe => VRPipe::StepOption->create(
                description   => 'path to your STAR executable',
                optional      => 1,
                default_value => 'STAR'
            ),
            star_sample_name_from_metadata => VRPipe::StepOption->create(
                description   => 'sample name from metadata to use in bam read groups',
                optional      => 1,
                default_value => 'sample'
            ),
            reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file'),
        };
    }
    
    method inputs_definition {
        return {
            fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => -1,
                description => 'fastq file(s) to be mapped',
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
                    optional       => ['mate', 'library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study']
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $star_exe      = $options->{star_exe};
            my $star_map_opts = $options->{star_map_options};
            my $star_gen_opts = $options->{star_genomeGenerate_options};
            my $star_sample   = $options->{star_sample_name_from_metadata};
            my $ref           = file($options->{reference_fasta});
            
            if ($star_map_opts =~ /runMode|genomeDir|readFilesIn|outSAMtype/) {
                $self->throw("star_map_options should not include the input fastq(s) and genome directory path");
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'STAR', version => 0, summary => 'STAR --genomeDir ' . $ref->dir . ' --readFilesIn fastqs ' . $star_map_opts));
            
            my %fastqs;
            foreach my $fq (@{ $self->inputs->{fastq_files} }) {
                my $meta = $fq->metadata;
                next if ($meta->{paired} == 0); #STAR can only handle the input of 2 paired-end fastqs
                my $lane = $meta->{lane};
                if ($meta->{paired} == 1) {
                    unshift @{ $fastqs{$lane} }, $fq;
                }
                else {
                    push @{ $fastqs{$lane} }, $fq;
                }
            }
            
            foreach my $lane (keys %fastqs) {
                my @fqs = map { $_->path } @{ $fastqs{$lane} };
                my $fq_meta = $self->common_metadata($fastqs{$lane});
                # add metadata and construct readgroup info
                my $bam_meta = {
                    lane   => $lane,
                    paired => 1,
                };
                my $rg_arg = 'ID:' . $lane;
                if (defined $fq_meta->{library}) {
                    my $lb = $fq_meta->{library};
                    $bam_meta->{library} = $lb;
                    $rg_arg .= ' LB:' . $self->command_line_safe_string($lb);
                }
                if (defined $fq_meta->{$star_sample}) {
                    my $sm = $fq_meta->{$star_sample};
                    $bam_meta->{sample} = $sm;
                    $rg_arg .= ' SM:' . $self->command_line_safe_string($sm);
                }
                if (defined $fq_meta->{insert_size}) {
                    my $pi = $fq_meta->{insert_size};
                    $bam_meta->{insert_size} = $pi;
                    $rg_arg .= ' PI:' . sprintf("%0.0f", $pi);
                }
                if (defined $fq_meta->{center_name}) {
                    my $cn = $fq_meta->{center_name};
                    $bam_meta->{center_name} = $cn;
                    $rg_arg .= ' CN:' . $self->command_line_safe_string($cn);
                }
                if (defined $fq_meta->{platform}) {
                    my $pl = $fq_meta->{platform};
                    $bam_meta->{platform} = $pl;
                    $rg_arg .= ' PL:' . $self->command_line_safe_string($pl);
                }
                if (defined $fq_meta->{study}) {
                    my $ds = $fq_meta->{study};
                    $bam_meta->{study} = $ds;
                    $rg_arg .= ' DS:' . $self->command_line_safe_string($ds);
                }
                if (defined $fq_meta->{analysis_group}) {
                    $bam_meta->{analysis_group} = $fq_meta->{analysis_group};
                }
                if (defined $fq_meta->{population}) {
                    $bam_meta->{population} = $fq_meta->{population};
                }
                
                my $bam_file = $self->output_file(output_key => 'bam_files', basename => "Aligned.out.bam", type => 'bam', metadata => $bam_meta);
                $self->output_file(output_key => 'tab_file', basename => 'SJ.out.tab', type => 'txt');
                foreach (qw(Log.progress.out Log.out Log.final.out)) {
                    $self->output_file(output_key => 'log_files', basename => $_, type => 'txt');
                }
                
                my ($cpus)         = $star_gen_opts =~ m/--runThreadN\s*(\d+)/;
                my ($SAindexBases) = $star_gen_opts =~ m/--genomeSAindexNbases\s*(\d+)/;
                my ($buffer)       = $star_gen_opts =~ m/--limitIObufferSize\s*(\d+)/;
                my $size           = -s $ref;
                my $memory         = int((10 * $size + 6 * 4**$SAindexBases + $buffer * $cpus) / 1024 / 1024); #MB
                my $req = $self->new_requirements(memory => $memory, time => 1200, $cpus ? (cpus => $cpus) : ());
                
                my $cmd = $star_exe . " --genomeDir " . $ref->dir . " --readFilesIn " . join(' ', @fqs) . " --outSAMattrRGline " . $rg_arg . " --outSAMtype BAM Unsorted " . $star_map_opts;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::star_map_fastq', 'map_and_check', [$cmd, $req, { output_files => [$bam_file] }]);
            }
        
        };
    }
    
    method outputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'mapped bam file per endedness and lane',
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
                    optional       => ['mate', 'library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study', 'bases']
                }
            ),
            tab_file => VRPipe::StepIODefinition->create(
                type            => 'txt',
                description     => 'tab-delimited file containing collapsed splice junctions',
                check_existence => 0
            ),
            log_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'log files'
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Maps the input fastq(s) with STAR";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method map_and_check (ClassName|Object $self: Str $cmd_line) {
        my $bam_file = VRPipe::File->get(path => file('Aligned.out.bam')->absolute);
        
        my @fq_paths = my ($fq_1_path, $fq_2_path) = $cmd_line =~ /--readFilesIn (\S+) (\S+)/;
        ($fq_1_path && $fq_2_path) || $self->throw("cmd_line [$cmd_line] had no --readFilesIn output specified");
        
        my $expected_reads = 0;
        foreach (@fq_paths) {
            my $fq_file = VRPipe::File->get(path => $_);
            my $meta = $fq_file->metadata;
            $expected_reads += $meta->{reads};
        }
        
        $bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $bam_file->update_stats_from_disc(retries => 3);
        
        my $actual_reads = $bam_file->num_records;
        $bam_file->add_metadata({ reads => $actual_reads });
        
        if ($cmd_line =~ /outSAMunmapped\s+Within/) {
            if ($actual_reads == $expected_reads) {
                return 1;
            }
            else {
                $bam_file->unlink;
                $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the bam file, yet there were $expected_reads reads in the fastq file(s)");
            }
        }
        else {
            return 1; #don't do any checks
        }
    }
    
    method command_line_safe_string (Str $str) {
        # truncate at the first space, convert non-word chars to underscores
        $str =~ s/^(\S+)\s.*/$1/;
        $str =~ s/[^\w\-]/_/g;
        return $str;
    }
}

1;
