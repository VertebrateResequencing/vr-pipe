
=head1 NAME

VRPipe::Steps::bam_merge_lane_splits - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Steps::bam_merge_lane_splits with VRPipe::StepRole {
    method options_definition {
        return { bam_merge_keep_single_paired_separate => VRPipe::StepOption->create(description   => 'when merging bam files, separately merges single ended bam files and paired-end bam files, resulting in 2 merged bam files',
                                                                                     optional      => 1,
                                                                                     default_value => 1),
                 samtools_exe => VRPipe::StepOption->create(description   => 'path to your samtools executable',
                                                            optional      => 1,
                                                            default_value => 'samtools') };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type        => 'bam',
                                                               max_files   => -1,
                                                               description => '1 or more bam files to merge',
                                                               metadata    => {
                                                                             lane           => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                             library        => 'library name',
                                                                             sample         => 'sample name',
                                                                             center_name    => 'center name',
                                                                             platform       => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                             study          => 'name of the study',
                                                                             insert_size    => 'expected (mean) insert size if paired',
                                                                             analysis_group => 'project analysis group',
                                                                             population     => 'sample population',
                                                                             bases          => 'total number of base pairs',
                                                                             reads          => 'total number of reads (sequences)',
                                                                             paired         => '0=unpaired reads were mapped; 1=paired reads were mapped',
                                                                             mapped_fastqs  => 'comma separated list of the fastq file(s) that were mapped',
                                                                             chunk          => 'mapped_fastq(s) are this chunk of original fastq(s)',
                                                                             optional       => ['library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study'] }),
                 dict_file => VRPipe::StepIODefinition->create(type        => 'txt',
                                                               description => 'a sequence dictionary file for your reference fasta') };
    }
    
    method body_sub {
        return sub {
            my $self      = shift;
            my $options   = $self->options;
            my $samtools  = $options->{samtools_exe};
            my $separate  = $options->{bam_merge_keep_single_paired_separate};
            my $dict_path = $self->inputs->{dict_file}->[0]->path;
            
            my ($lane, %bams, %metas);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $this_path = $bam->path;
                
                my $meta      = $bam->metadata;
                my $this_lane = $meta->{lane};
                $lane ||= $this_lane;
                $self->throw("Not all the input bams were for the same lane") if $this_lane ne $lane;
                
                my $paired = $meta->{paired};
                push(@{ $bams{$paired} }, $this_path);
                
                unless (defined $metas{$paired}) {
                    $metas{$paired} = {};
                    foreach my $key (qw(lane library sample center_name platform study insert_size analysis_group population)) {
                        if (defined $meta->{$key}) {
                            $metas{$paired}->{$key} = $meta->{$key};
                        }
                    }
                    
                    my @fqs = split(',', $meta->{mapped_fastqs});
                    my @these_parents;
                    foreach my $fq (@fqs) {
                        my $fq_file = VRPipe::File->get(path => $fq, auto_resolve => 1);
                        my $parent = $fq_file->metadata->{source_fastq} || $self->throw("no source_fastq for one of the mapped fastqs of $this_path - was it really mapped from a fastq_split result?");
                        push(@these_parents, $parent);
                    }
                    $metas{$paired}->{mapped_fastqs} = join(',', @these_parents);
                }
                $metas{$paired}->{reads} += $meta->{reads};
                $metas{$paired}->{bases} += $meta->{bases};
                $metas{$paired}->{paired} = $paired;
            }
            
            unless ($separate) {
                if (keys %bams == 2) {
                    push(@{ $bams{2} }, @{ $bams{0} }, @{ $bams{1} });
                    delete $bams{0};
                    delete $bams{1};
                    
                    $metas{2} = $metas{1};
                    $metas{2}->{reads} += $metas{0}->{reads};
                    $metas{2}->{bases} += $metas{0}->{bases};
                    $metas{2}->{mapped_fastqs} = join(',', $metas{0}->{mapped_fastqs}, $metas{1}->{mapped_fastqs});
                    $metas{2}->{paired} = 2;
                }
            }
            
            my $req = $self->new_requirements(memory => 1000, time => 1);
            my $step_state = $self->step_state->id;
            while (my ($paired, $in_bams) = each %bams) {
                my $basename = $lane;
                if ($separate) {
                    $basename .= $paired == 0 ? '.se' : '.pe';
                }
                $basename .= '.bam';
                my $merge_file = $self->output_file(output_key => 'merged_lane_bams',
                                                    basename   => $basename,
                                                    type       => 'bam',
                                                    metadata   => $metas{$paired});
                my $merge_path = $merge_file->path;
                
                $self->output_file(basename  => $basename . '.header',
                                   type      => 'txt',
                                   temporary => 1);
                
                if (@$in_bams == 1) {
                    my $sam_file = $basename;
                    $sam_file =~ s/\.bam$/.sam/;
                    $self->output_file(basename  => $sam_file,
                                       type      => 'txt',
                                       temporary => 1);
                }
                
                my $this_cmd = "use VRPipe::Steps::bam_merge_lane_splits; VRPipe::Steps::bam_merge_lane_splits->merge_and_check(samtools => q[$samtools], dict => q[$dict_path], output => q[$merge_path], step_state => $step_state, bams => [qw(@$in_bams)]);";
                $self->dispatch_vrpipecode($this_cmd, $req); # deliberately do not include {output_files => [$merge_file]} so that any temp files we made will get their stats updated prior to auto-deletion
            }
        };
    }
    
    method outputs_definition {
        return { merged_lane_bams => VRPipe::StepIODefinition->create(type        => 'bam',
                                                                      max_files   => 2,
                                                                      description => 'a merged bam file for each library layout (single ended vs paired)',
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
                                                                                    paired         => '0=unpaired reads were mapped; 1=paired reads were mapped; 2=mixture of paired and unpaired reads were mapped',
                                                                                    mapped_fastqs  => 'comma separated list of the fastq file(s) that were mapped',
                                                                                    optional       => ['library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study'] }) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merges multiple bam files for the same lane (mapped from splits of the same fastq pair) into a single bam file (per library layout). Also ensures the header has complete sequence information, a good RG line, and chained PG lines"; #*** , and that all records have an RG tag
    }
    
    method max_simultaneous {
        return 0;                                                                                                                                                                                                                                      # meaning unlimited
    }
    
    method merge_and_check (ClassName|Object $self: Str|File :$samtools!, Str|File :$dict!, Str|File :$output!, Persistent :$step_state!, ArrayRef[Str|File] :$bams!) {
        # make a nice sam header
        my $header_file = VRPipe::File->get(path => $output . '.header');
        my $header_path = $header_file->path;
        my $hfh         = $header_file->openw;
        
        my $header_lines = 0;
        print $hfh "\@HD\tVN:1.0\tSO:coordinate\n";
        $header_lines++;
        
        # copy over the SQ lines from the dict file
        my $dict_file = VRPipe::File->get(path => $dict);
        my $dfh = $dict_file->openr;
        while (<$dfh>) {
            next unless /^\@SQ/;
            print $hfh $_;
            $header_lines++;
        }
        
        # construct the RG line from the bam metadata
        my $merged_bam_file = VRPipe::File->get(path => $output);
        my $meta = $merged_bam_file->metadata;
        print $hfh "\@RG\tID:", $meta->{lane};
        if (defined $meta->{library}) {
            print $hfh "\tLB:", $meta->{library};
        }
        if (defined $meta->{sample}) {
            print $hfh "\tSM:", $meta->{sample};
        }
        if (defined $meta->{insert_size}) {
            print $hfh "\tPI:", $meta->{insert_size};
        }
        if (defined $meta->{center_name}) {
            print $hfh "\tCN:", $meta->{center_name};
        }
        if (defined $meta->{platform}) {
            print $hfh "\tPL:", $meta->{platform};
        }
        if (defined $meta->{study}) {
            print $hfh "\tDS:", $meta->{study};
        }
        print $hfh "\n";
        $header_lines++;
        
        # construct a chain of PG lines for the header by looking at previous
        # steps in our pipeline
        my $this_step_state = VRPipe::StepState->get(id => $step_state);
        my $pipelinesetup   = $this_step_state->pipelinesetup;
        my $dataelement     = $this_step_state->dataelement;
        my $stepmember      = $this_step_state->stepmember;
        my $this_stepm_id   = $stepmember->id;
        my $pipeline        = $stepmember->pipeline;
        my $pp;
        my %step_name_counts;
        
        foreach my $stepm ($pipeline->step_members) {
            last if $stepm->id == $this_stepm_id;
            
            my $cmd_summary = VRPipe::StepState->get(pipelinesetup => $pipelinesetup, stepmember => $stepm, dataelement => $dataelement)->cmd_summary || next;
            
            my $step_name = $stepm->step->name;
            my $snc       = ++$step_name_counts{$step_name};
            if ($snc > 1) {
                $step_name .= '.' . $snc;
            }
            print $hfh "\@PG\tID:$step_name\tPN:", $cmd_summary->exe, "\t";
            if ($pp) {
                print $hfh "PP:$pp\t";
            }
            print $hfh "VN:", $cmd_summary->version, "\tCL:", $cmd_summary->summary, "\n";
            $header_lines++;
            
            $pp = $step_name;
        }
        $header_file->close;
        
        my $cmd_line;
        if (@$bams > 1) {
            $cmd_line = "$samtools merge -h $header_path $output @$bams";
        }
        else {
            my $sam_path = $output;
            $sam_path =~ s/\.bam$//;
            $sam_path .= '.sam';
            my $bam_path = $bams->[0];
            $cmd_line = "cat $header_path > $sam_path; $samtools view $bam_path >> $sam_path; $samtools view -bS $sam_path > $output";
        }
        $merged_bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $expected_lines = $meta->{reads} + $header_lines;
        $merged_bam_file->update_stats_from_disc(retries => 3);
        my $actual_lines = $merged_bam_file->lines;
        
        if ($actual_lines == $expected_lines) {
            return 1;
        }
        else {
            $merged_bam_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_lines lines were generated in the merged bam file, yet there were $expected_lines records in the input bam files and header");
        }
    }
}

1;
