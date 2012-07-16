=head1 NAME

VRPipe::Steps::stampy_map_fastq - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Steps::stampy_map_fastq with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file to map against'),
                 stampy_map_options => VRPipe::StepOption->create(description => 'options for stampy mapping, excluding the output, input fastq(s), reference, --bwa and --bwaoptions options',
                                                               optional => 1),
                 stampy_bwa_options => VRPipe::StepOption->create(description => 'to use bwa for premapping, supply the string you would give to --bwaoptions, excluding the reference',
                                                               optional => 1),
                 stampy_substitution_rate_from_metadata => VRPipe::StepOption->create(description => q[set stampy's --substitutionrate option based on the value of the substitution_rate metadata key on the input fastqs, if present],
                                                                                   optional => 1,
                                                                                   default_value => 1),
                 stampy_exe => VRPipe::StepOption->create(description => 'path to your stampy.py executable',
                                                       optional => 1,
                                                       default_value => 'stampy.py'),
                 bwa_exe => VRPipe::StepOption->create(description => 'path to your bwa executable',
                                                    optional => 1,
                                                    default_value => 'bwa') };
    }
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->create(type => 'fq',
                                                              max_files => -1,
                                                              description => 'fastq file(s) to be mapped',
                                                              metadata => {lane => 'lane name (a unique identifer for this sequencing run)',
                                                                           library => 'library name',
                                                                           sample => 'sample name',
                                                                           center_name => 'center name',
                                                                           platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                           study => 'name of the study, put in the DS field of the RG header line',
                                                                           insert_size => 'expected (mean) insert size if paired',
                                                                           analysis_group => 'project analysis group',
                                                                           population => 'sample population',
                                                                           bases => 'total number of base pairs',
                                                                           reads => 'total number of reads (sequences)',
                                                                           paired => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                                                                           mate => 'if paired, the path to the fastq that is our mate',
                                                                           chunk => 'if the fastq file was produced by fastq_split Step, the chunk number',
                                                                           substitution_rate => 'the substitution rate of this sample data compared to the reference',
                                                                           optional => ['mate', 'chunk', 'library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study', 'substitution_rate']}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $stampy_exe = $options->{stampy_exe};
            my $stampy_opts = $options->{stampy_map_options};
            if ($stampy_opts =~ /$ref|-h|-g|--bwa=|bwaoptions|-M|-o|-f/) {
                $self->throw("stampy_map_options should not include the output, input fastq(s), reference, --bwa nor --bwaoptions options");
            }
            
            my $bwa_opts = $options->{stampy_bwa_options};
            if ($bwa_opts) {
                if ($bwa_opts =~ /$ref/) {
                    $self->throw("stampy_bwa_options should not include the reference");
                }
                $stampy_opts .= ' --bwa='.$options->{bwa_exe}.' --bwaoptions={'.$bwa_opts.' '.$ref.'}'; # {} instead of "" because of quoting in shell issues
            }
            
            my $req = $self->new_requirements(memory => 5900, time => 2);
            my $cmd = $stampy_exe.' '.$stampy_opts." -g $ref -h $ref ";
            
            # we can handle the input of both a single-ended fastq, and 2
            # paired-end fastqs, and groups of these (up to) 3 files for
            # multiple different lanes
            my @fq_files = @{$self->inputs->{fastq_files}};
            my %fqs_by_lane;
            my %fqs_by_path;
            foreach my $fq (@fq_files) {
                my $fq_meta = $fq->metadata;
                my $paired = $fq_meta->{paired};
                my $path = $fq->resolve->path->stringify; # since we're storing this path in output metadata, we want the real path of the fq, not that of a symlink
                my $lane = $fq_meta->{lane};
                my $chunk = $fq_meta->{chunk} || 0;
                $fqs_by_path{$path} = [$lane, $chunk, $paired, $fq_meta];
                $fqs_by_lane{$lane}->{$chunk}->{$paired} = $path;
            }
            
            my %srates;
            while (my ($lane, $chunks) = each %fqs_by_lane) {
                while (my ($chunk, $ends) = each %$chunks) {
                    while (my ($paired, $path) = each %$ends) {
                        next if $paired == 2;
                        
                        my @fqs = ($path);
                        my $fq_meta = $fqs_by_path{$fqs[0]}->[3];
                        my $reads = $fq_meta->{reads};
                        my $bases = $fq_meta->{bases};
                        
                        my $this_cmd;
                        if ($paired == 0) {
                            $this_cmd = $cmd."-M $fqs[0]";
                        }
                        else {
                            my $mate_path = $ends->{2};
                            push(@fqs, $mate_path);
                            
                            $this_cmd = $cmd."-M $fqs[0],$fqs[1]";
                            
                            $reads += $fqs_by_path{$fqs[1]}->[3]->{reads};
                            $bases += $fqs_by_path{$fqs[1]}->[3]->{bases};
                        }
                        
                        my $srate = $fq_meta->{substitution_rate};
                        unless (defined $srate) {
                            # also check the source fastq or bam, in case this
                            # is a fastq split or was made from a bam
                            SOURCE: foreach my $key (qw(source_bam source_fastq)) {
                                my @source_paths = split(',', $fq_meta->{$key} || next);
                                foreach my $path (@source_paths) {
                                    my $parent_srate = VRPipe::File->get(path => $path)->metadata->{substitution_rate};
                                    if (defined $parent_srate) {
                                        $srate = $parent_srate;
                                        last SOURCE;
                                    }
                                }
                            }
                        }
                        if ($options->{stampy_substitution_rate_from_metadata} && defined $srate) {
                            if ($srate == 0) {
                                # stampy doesn't like 0, and the default is 0.001
                                $srate = '0.00001';
                            }
                            $this_cmd =~ s/--substitutionrate[= ]+\S+//;
                            $this_cmd .= ' --substitutionrate='.$srate;
                            $srates{$srate} = 1;
                        }
                        
                        # add metadata and construct readgroup info
                        my $rg_arg = 'ID:'.$lane;
                        my $sam_meta = {lane => $lane,
                                        bases => $bases,
                                        reads => $reads,
                                        paired => $paired,
                                        mapped_fastqs => join(',', @fqs),
                                        $chunk ? (chunk => $chunk) : ()};
                        if (defined $fq_meta->{library}) {
                            my $lb = $fq_meta->{library};
                            $sam_meta->{library} = $lb;
                            $rg_arg .= ',LB:'.$self->command_line_safe_string($lb);
                        }
                        if (defined $fq_meta->{sample}) {
                            my $sm = $fq_meta->{sample};
                            $sam_meta->{sample} = $sm;
                            $rg_arg .= ',SM:'.$self->command_line_safe_string($sm);
                        }
                        if (defined $fq_meta->{insert_size}) {
                            my $pi = $fq_meta->{insert_size};
                            $sam_meta->{insert_size} = $pi;
                            $rg_arg .= ',PI:'.sprintf("%0.0f", $pi); # stampy enforces an int for PL
                        }
                        if (defined $fq_meta->{center_name}) {
                            my $cn = $fq_meta->{center_name};
                            $sam_meta->{center_name} = $cn;
                            $rg_arg .= ',CN:'.$self->command_line_safe_string($cn);
                        }
                        if (defined $fq_meta->{platform}) {
                            my $pl = $fq_meta->{platform};
                            $sam_meta->{platform} = $pl;
                            $rg_arg .= ',PL:'.$self->command_line_safe_string($pl);
                        }
                        if (defined $fq_meta->{study}) {
                            my $ds = $fq_meta->{study};
                            $sam_meta->{study} = $ds;
                            $rg_arg .= ',DS:'.$self->command_line_safe_string($ds);
                        }
                        if (defined $fq_meta->{analysis_group}) {
                            $sam_meta->{analysis_group} = $fq_meta->{analysis_group};
                        }
                        if (defined $fq_meta->{population}) {
                            $sam_meta->{population} = $fq_meta->{population};
                        }
                        
                        my $ended = $paired ? 'pe' : 'se';
                        my $sam_file = $self->output_file(output_key => 'stampy_sam_files',
                                                          basename => $chunk ? "$lane.$ended.$chunk.sam" : "$lane.$ended.sam",
                                                          type => 'txt',
                                                          metadata => $sam_meta);
                        
                        $this_cmd .= " --readgroup=$rg_arg -o ".$sam_file->path;
                        $self->dispatch_wrapped_cmd('VRPipe::Steps::stampy_map_fastq', 'map_and_check', [$this_cmd, $req, {output_files => [$sam_file]}]);
                    }
                }
            }
            
            my @srates = keys %srates;
            if (@srates == 1) {
                $stampy_opts =~ s/--substitutionrate[= ]+\S+//;
                $stampy_opts .= ' ' if $stampy_opts;
                $stampy_opts .= '--substitutionrate='.$srates[0];
            }
            $stampy_opts =~ s/$ref/\$ref/g;
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'stampy', version => VRPipe::StepCmdSummary->determine_version($stampy_exe, '^stampy v(\S+)'), summary => 'stampy.py '.$stampy_opts.' -g $ref.fa -h $ref.fa -o $out.sam -M $fastq(s)'));
        };
    }
    method outputs_definition {
        return { stampy_sam_files => VRPipe::StepIODefinition->create(type => 'txt',
                                                                   max_files => -1,
                                                                   description => 'mapped sam file per endedness and lane',
                                                                   metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                                library => 'library name',
                                                                                sample => 'sample name',
                                                                                center_name => 'center name',
                                                                                platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                                study => 'name of the study, put in the DS field of the RG header line',
                                                                                insert_size => 'expected (mean) insert size if paired',
                                                                                analysis_group => 'project analysis group',
                                                                                population => 'sample population',
                                                                                bases => 'total number of base pairs',
                                                                                reads => 'total number of reads (sequences)',
                                                                                paired => '0=unpaired reads were mapped; 1=paired reads were mapped',
                                                                                mapped_fastqs => 'comma separated list of the fastq file(s) that were mapped',
                                                                                chunk => 'if this was mapped with fastqs that were chunks of an original fastq, this tells you which chunk',
                                                                                optional => ['chunk', 'library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study']}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Maps the input fastq(s) with stampy to the reference";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method map_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($sam_path) = $cmd_line =~ /-o (\S+)/;
        $sam_path || $self->throw("cmd_line [$cmd_line] had no -o output specified");
        $cmd_line =~ s/\{/"/;
        $cmd_line =~ s/\}/"/;
        
        my $sam_file = VRPipe::File->get(path => $sam_path);
        my $expected_reads = $sam_file->metadata->{reads};
        
        $sam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $sam_file->update_stats_from_disc(retries => 3);
        my $lines = $sam_file->lines;
        
        if ($lines > $expected_reads) {
            return 1;
        }
        else {
            $sam_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $lines lines were generated in the sam file, yet there were $expected_reads reads in the fastq file(s)");
        }
    }
    
    method command_line_safe_string (Str $str) {
        # truncate at the first space, convert non-word chars to underscores
        $str =~ s/^(\S+)\s.*/$1/;
        $str =~  s/[^\w\-]/_/g;
        return $str;
    }
}

1;
