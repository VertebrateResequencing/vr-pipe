=head1 NAME

VRPipe::Steps::bam_merge - a step

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

class VRPipe::Steps::bam_merge extends VRPipe::Steps::picard {
    around options_definition {
        return { %{$self->$orig},
                 merge_sam_files_options => VRPipe::StepOption->create(description => 'options for picard MergeSamFiles', optional => 1, default_value => 'VALIDATION_STRINGENCY=SILENT'),
                 bam_merge_keep_single_paired_separate => VRPipe::StepOption->create(description => 'when merging bam files, separately merges single ended bam files and paired-end bam files, resulting in 2 merged bam files',
                                                                                  optional => 1,
                                                                                  default_value => 1),
                 bam_merge_maximum_files => VRPipe::StepOption->create(description => 'set the maximum number of files to merge in one go; eg. if there are 100 bams to merge and you set this to 40, 3 output bams will result, being merges of 40, 40 and 20 of the input bams (defaults to infinity)',
                                                                    optional => 1) };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', 
                                                            max_files => -1, 
                                                            description => '1 or more bam files',
                                                            metadata => {lane => 'lane name, or comma separated list of lane names if merged',
                                                                         library => 'library name, or comma separated list of library names if merged',
                                                                         sample => 'sample name, or comma separated list of sample names if merged',
                                                                         center_name => 'center name, or comma separated list of center names if merged',
                                                                         platform => 'sequencing platform, or comma separated list of platform names if merged',
                                                                         study => 'name of the study, or comma separated list of study names if merged',
                                                                         bases => 'total number of base pairs',
                                                                         reads => 'total number of reads (sequences)',
                                                                         paired => '0=unpaired reads were mapped; 1=paired reads were mapped',
                                                                         optional => ['lane', 'library', 'sample', 'center_name', 'platform', 'study', 'bases']}
                                                            ),
                };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            my $merge_jar = $self->jar('MergeSamFiles.jar');
            
            my $opts = $options->{merge_sam_files_options};
            if ($opts =~ /MergeSamFiles/) {
                $self->throw("merge_sam_files_options should not include the MergeSamFiles task command");
            }
            
            my %paired_groups;
            my $paired = 0;
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $meta = $bam->metadata;
                $paired ||= $meta->{paired};
                my $ended = $meta->{paired} ? 'pe' : 'se';
                push @{$paired_groups{$ended}}, $bam;
            }
            
            # merge single and paired end bams unless option says not to
            unless ($options->{bam_merge_keep_single_paired_separate}) {
                if ($paired && exists $paired_groups{se}) {
                    push @{$paired_groups{pe}}, @{delete $paired_groups{se}};
                }
            }
            
            # break up the paired groups according to bam_merge_maximum_files
            my $max_bams_per_group = $options->{bam_merge_maximum_files};
            my %merge_groups;
            while (my ($ended, $bam_files) = each %paired_groups) {
                if ($max_bams_per_group) {
                    my @these_bams;
                    my $chunk = 1;
                    foreach my $bam (@$bam_files) {
                        push(@these_bams, $bam);
                        if (@these_bams == $max_bams_per_group) {
                            $merge_groups{$ended}->{$chunk++} = [@these_bams];
                            @these_bams = ();
                        }
                    }
                    if (@these_bams) {
                        $merge_groups{$ended}->{$chunk} = [@these_bams];
                    }
                }
                else {
                    $merge_groups{$ended}->{1} = $bam_files;
                }
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'picard', 
                                   version => $self->picard_version(),
                                   summary => 'java $jvm_args -jar MergeSamFiles.jar INPUT=$bam_file(s) OUTPUT=$merged_bam '.$opts));
            
            my $req = $self->new_requirements(memory => 1000, time => 1);
            while (my ($ended, $chunks) = each %merge_groups) {
                while (my ($chunk, $bam_files) = each %$chunks) {
                    my $merge_file = $self->output_file(output_key => 'merged_bam_files',
                                                        basename => "$ended.$chunk.bam",
                                                        type => 'bam',
                                                        metadata => $self->common_metadata($bam_files)
                                                        );
                    
                    my $temp_dir = $options->{tmp_dir} || $merge_file->dir;
                    my $jvm_args = $self->jvm_args($req->memory, $temp_dir);
                    
                    my $in_bams = join ' INPUT=', map { $_->path } @$bam_files;
                    $in_bams = ' INPUT='.$in_bams;
                    my $this_cmd = $self->java_exe.qq[ $jvm_args -jar $merge_jar$in_bams OUTPUT=].$merge_file->path.qq[ $opts];
                    $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_merge', 'merge_and_check', [$this_cmd, $req, {output_files => [$merge_file]}]);
                }
            }
        };
    }
    method outputs_definition {
        return { merged_bam_files => VRPipe::StepIODefinition->create(type => 'bam', 
                                                                   max_files => -1, 
                                                                   description => '1 or more bam merged bam files',
                                                                   metadata => {lane => 'lane name, or comma separated list of lane names if merged',
                                                                                library => 'library name, or comma separated list of library names if merged',
                                                                                sample => 'sample name, or comma separated list of sample names if merged',
                                                                                center_name => 'center name, or comma separated list of center names if merged',
                                                                                platform => 'sequencing platform, or comma separated list of platform names if merged',
                                                                                study => 'name of the study, or comma separated list of study names if merged',
                                                                                bases => 'total number of base pairs',
                                                                                reads => 'total number of reads (sequences)',
                                                                                paired => '0=unpaired reads were mapped; 1=paired reads were mapped',
                                                                                optional => ['lane', 'library', 'sample', 'center_name', 'platform', 'study', 'bases']}
                                                                    ),
               };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Merges bam files using Picard MergeSamFiles";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    method merge_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($out_path) = $cmd_line =~ /OUTPUT=(\S+)/;
        my @in_paths = $cmd_line =~ /INPUT=(\S+)/g;
        @in_paths || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $out_file = VRPipe::File->get(path => $out_path);
        my @in_files;
        foreach my $in_path (@in_paths) {
            push @in_files, VRPipe::File->get(path => $in_path);
        }
        
        if (scalar @in_paths == 1) {
            $in_files[0]->copy($out_file);
            return;
        }
        else {
            $out_file->disconnect;
            system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        }
        
        $out_file->update_stats_from_disc(retries => 3);
        
        my $reads = 0;
        my $bases = 0;
        my $paired = 0;
        foreach my $in_file (@in_files) {
            my $meta = $in_file->metadata;
            $reads += $meta->{reads};
            $bases += $meta->{bases} if $meta->{bases};
            $paired ||= $meta->{paired};
        }
        my $actual_reads = $out_file->num_records;
        
        if ($actual_reads == $reads) {
            my %new_meta;
            $new_meta{reads} = $actual_reads;
            $new_meta{bases} = $bases if $bases;
            $new_meta{paired} = $paired;
            $out_file->add_metadata(\%new_meta);
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output bam file, yet there were $reads reads in the original bam files");
        }
    }
}

1;
