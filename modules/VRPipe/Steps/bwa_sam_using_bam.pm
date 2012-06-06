=head1 NAME

VRPipe::Steps::bwa_sam_using_bam - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Steps::bwa_sam_using_bam with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file the sai files were aligned with'),
                 bwa_samse_options => VRPipe::StepOption->get(description => 'options to bwa samse, excluding the input sai, bam, reference, -r and -f',
                                                              optional => 1,
                                                              default_value => ''),
                 bwa_sampe_options => VRPipe::StepOption->get(description => 'options to bwa sampe, excluding the input sai, bams, reference, -r and -f; defaults to bwa sampe with -a set as appropriate for the bam insert_size',
                                                              optional => 1),
                 bwa_exe => VRPipe::StepOption->get(description => 'path to your bwa executable',
                                                    optional => 1,
                                                    default_value => 'bwa') };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam',
                                                            max_files => -1,
                                                            description => 'bam files',
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
                                                                         forward_reads => 'number of forward reads',
                                                                         reverse_reads => 'number of reverse reads',
                                                                         paired => '0=single ended reads only; 1=paired end reads present',
                                                                         optional => ['library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study']}),
                 sai_files => VRPipe::StepIODefinition->get(type => 'bin',
                                                            max_files => -1,
                                                            description => 'sai files, as produced by bwa aln',
                                                            metadata => {source_bam => 'the bam file that was input to bwa aln to generate this sai file',
                                                                         paired => 'what kind of reads were used from the source_bam: 0=single end reads; 1=forward reads; 2=reverse reads',
                                                                         reference_fasta => 'the reference fasta that reads were aligned to'}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my %cmds;
            my $bwa_exe = $options->{bwa_exe};
            foreach my $ended ('se', 'pe') {
                my $opts = $options->{"bwa_sam${ended}_options"} || '';
                if ($opts =~ /$ref|-r|-f|sam$ended/) {
                    $self->throw("bwa_sam${ended}_cmd should not include the ref or -r or -f option, or the sam$ended sub command");
                }
                
                my $cmd = $bwa_exe." sam$ended ".$opts;
                $cmd =~ s/\s+$//;
                $cmds{$ended} = $cmd;
            }
            
            my $req = $self->new_requirements(memory => 4900, time => 4);
            
            my @bam_files = @{$self->inputs->{bam_files}};
            my @sai_files = @{$self->inputs->{sai_files}};
            
            my %bams_by_path;
            foreach my $bam (@bam_files) {
                my $bam_meta = $bam->metadata;
                my $path = $bam->path->stringify;
                my $lane = $bam_meta->{lane};
                $bams_by_path{$path} = [$lane, $bam_meta];
            }
            
            my %by_lane;
            foreach my $sai (@sai_files) {
                my $sai_meta = $sai->metadata;
                my $path = $sai->path->stringify;
                my $source_bam = $sai_meta->{source_bam};
                my $read_type = $sai_meta->{paired};
                my $paired = $read_type == 0 ? 'se' : 'pe';
                my $ref = $bams_by_path{$source_bam};
                $ref || $self->throw("got a sai file $path with source_bam $source_bam, but that bam was not an input fastq to bwa_sam_using_bam Step");
                push(@{$by_lane{$ref->[0]}->{$paired}->{sais}}, $path);
                $by_lane{$ref->[0]}->{$paired}->{bam} = $source_bam;
            }
            
            my $summary_cmd;
            while (my ($lane, $ends) = each %by_lane) {
                while (my ($paired, $paths) = each %$ends) {
                    my $bam = $paths->{bam};
                    my @sais = @{$paths->{sais}};
                    my $bam_meta = $bams_by_path{$bam}->[1];
                    
                    # don't create an se sam if it will have no reads in it
                    my $paired_reads = ($bam_meta->{forward_reads} || 0) + ($bam_meta->{reverse_reads} || 0);
                    if ($paired eq 'se' && $bam_meta->{paired}) {
                        next if $paired_reads == $bam_meta->{reads};
                    }
                    
                    # don't create a pe sam if it's all single end reads
                    if ($paired eq 'pe' && ($paired_reads == 0 || ! $bam_meta->{paired})) {
                        next;
                    }
                    
                    my $reads = $bam_meta->{reads}; #*** what happens when a bam had mixed paired and unpaired reads?
                    my $bases = $bam_meta->{bases};
                    
                    my $this_cmd;
                    my @bams = ($bam);
                    if ($paired eq 'se') {
                        $this_cmd = $cmds{se};
                        $summary_cmd ||= $this_cmd;
                    }
                    else {
                        $this_cmd = $cmds{pe};
                        push(@bams, $bam);
                        
                        unless ($this_cmd =~ /-a/) {
                            my $insert_size = $bam_meta->{insert_size} || 500;
                            my $max = $insert_size * 3;
                            $this_cmd .= " -a $max";
                        }
                        $summary_cmd = $this_cmd;
                    }
                    
                    # add metadata and construct RG line
                    my $rg_line = '@RG\tID:'.$lane;
                    my $sam_meta = {lane => $lane,
                                    bases => $bases,
                                    reads => $reads,
                                    paired => $paired eq 'se' ? 0 : 1,
                                    remapped_bam => $bam};
                    if (defined $bam_meta->{library}) {
                        my $lb = $bam_meta->{library};
                        $sam_meta->{library} = $lb;
                        $rg_line .= '\tLB:'.$lb;
                    }
                    if (defined $bam_meta->{sample}) {
                        my $sm = $bam_meta->{sample};
                        $sam_meta->{sample} = $sm;
                        $rg_line .= '\tSM:'.$sm;
                    }
                    if (defined $bam_meta->{insert_size}) {
                        my $pi = $bam_meta->{insert_size};
                        $sam_meta->{insert_size} = $pi;
                        $rg_line .= '\tPI:'.$pi;
                    }
                    if (defined $bam_meta->{center_name}) {
                        my $cn = $bam_meta->{center_name};
                        $sam_meta->{center_name} = $cn;
                        $rg_line .= '\tCN:'.$cn;
                    }
                    if (defined $bam_meta->{platform}) {
                        my $pl = $bam_meta->{platform};
                        $sam_meta->{platform} = $pl;
                        $rg_line .= '\tPL:'.$pl;
                    }
                    if (defined $bam_meta->{study}) {
                        my $ds = $bam_meta->{study};
                        $sam_meta->{study} = $ds;
                        $rg_line .= '\tDS:'.$ds;
                    }
                    if (defined $bam_meta->{analysis_group}) {
                        $sam_meta->{analysis_group} = $bam_meta->{analysis_group};
                    }
                    if (defined $bam_meta->{population}) {
                        $sam_meta->{population} = $bam_meta->{population};
                    }
                    
                    my $lane_base = $lane;
                    $lane_base =~ s/\W/_/g;
                    my $sam_file = $self->output_file(output_key => 'bwa_sam_files',
                                                      basename => "$lane_base.$paired.sam",
                                                      type => 'txt',
                                                      metadata => $sam_meta);
                    
                    $this_cmd .= " -r '$rg_line' -f ".$sam_file->path." $ref @sais @bams";
                    $self->dispatch_wrapped_cmd('VRPipe::Steps::bwa_sam_using_bam', 'sam_and_check', [$this_cmd, $req, {output_files => [$sam_file]}]);
                }
            }
            
            $summary_cmd =~ s/^\S+/bwa/;
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'bwa', version => VRPipe::StepCmdSummary->determine_version($bwa_exe, '^Version: (.+)$'), summary => $summary_cmd.' -r $rg_line -f $sam_file $reference_fasta $sai_file(s) $bam_file(s)'));
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
                                                                             analysis_group => 'project analysis group',
                                                                             population => 'sample population',
                                                                             bases => 'total number of base pairs',
                                                                             reads => 'total number of reads (sequences)',
                                                                             paired => '0=unpaired reads were mapped; 1=paired reads were mapped',
                                                                             remapped_bam => 'the bam file that was (re)mapped',
                                                                             optional => ['library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study']}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Produces sam files with bwa samse/sampe for each input bams&sais";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method sam_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($sam_path) = $cmd_line =~ /-f (\S+)/;
        $sam_path || $self->throw("cmd_line [$cmd_line] had no -f output specified");
        
        my $sam_file = VRPipe::File->get(path => $sam_path);
        
        $sam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]"); #*** want to exit with the exit code of bwa failing...
        
        my $expected_reads = $sam_file->metadata->{reads};
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
}

1;
