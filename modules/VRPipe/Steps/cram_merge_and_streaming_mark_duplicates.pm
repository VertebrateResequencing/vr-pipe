
=head1 NAME

VRPipe::Steps::cram_merge_and_streaming_mark_duplicates - a step

=head1 DESCRIPTION

Runs the samtools merge and biobambam's bamstreamingmarkduplicates to mark
duplicates.

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::cram_merge_and_streaming_mark_duplicates with VRPipe::StepRole {
    method options_definition {
        return {
            samtools_exe                    => VRPipe::StepOption->create(description => 'path to samtools executable',                                                               optional => 1, default_value => 'samtools'),
            bamstreamingmarkduplicates_exe  => VRPipe::StepOption->create(description => 'path to bamstreamingmarkduplicates executable',                                             optional => 1, default_value => 'bamstreamingmarkduplicates'),
            bamstreamingmarkduplicates_opts => VRPipe::StepOption->create(description => 'bamstreamingmarkduplicates options (excluding arguments that set input/output file names)', optional => 1, default_value => 'resetdupflag=1'),
            tmp_dir                         => VRPipe::StepOption->create(description => 'location for tmp directories; defaults to working directory',                               optional => 1),
        };
    }
    
    method inputs_definition {
        return {
            cram_files => VRPipe::StepIODefinition->create(
                type        => 'cram',
                max_files   => -1,
                description => '1 or more coordinate sorted CRAM files',
                metadata    => {
                    bases    => 'total number of base pairs',
                    reads    => 'total number of reads (sequences)',
                    paired   => '0=single ended reads only; 1=paired end reads present',
                    optional => ['bases', 'reads', 'paired']
                }
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self                            = shift;
            my $options                         = $self->options;
            my $samtools                        = $options->{samtools_exe};
            my $bamstreamingmarkduplicates_exe  = $options->{bamstreamingmarkduplicates_exe};
            my $bamstreamingmarkduplicates_opts = $options->{bamstreamingmarkduplicates_opts};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bamstreamingmarkduplicates',
                    version => VRPipe::StepCmdSummary->determine_version($bamstreamingmarkduplicates_exe . ' --version', '^This is biobambam version (.+)\.$'),
                    summary => "bamstreamingmarkduplicates $bamstreamingmarkduplicates_opts I=\$cram_file(s) O=\$output_file"
                )
            );
            
            # my $req = $self->new_requirements(memory => 3000, time => 1, cpus => $bamstreamingmarkduplicates_threads);
            my $req = $self->new_requirements(memory => 3000, time => 1);
            
            my $inputs           = $self->inputs->{cram_files};
            my $merged_metadata  = $self->common_metadata($inputs);
            my @input_paths      = map { $_->path } @$inputs;
            my $file_list        = VRPipe::FileList->create(files => $inputs)->id;
            my $basename         = $self->step_state->dataelement->id;
            my $markdup_bam_file = $self->output_file(
                output_key => 'markdup_bam_files',
                basename   => $basename . '.bam',
                type       => 'bam',
                metadata   => $merged_metadata
            );
            my $markdup_bai_file = $self->output_file(
                output_key => 'markdup_bai_files',
                basename   => $basename . '.bam.bai',
                type       => 'bai',
                metadata   => $merged_metadata
            );
            my $markdup_metrics_file = $self->output_file(
                output_key => 'markdup_metrics_files',
                basename   => $basename . '.metrics',
                type       => 'txt',
                metadata   => $merged_metadata
            );
            my $tmpdir = $options->{tmp_dir} || $markdup_bam_file->dir;
            my $preprocess_cmd = scalar @input_paths > 1 ? qq[$samtools merge -u - @input_paths] : qq[$samtools view -u @input_paths];
            my $markdup_cmd    = "$bamstreamingmarkduplicates_exe $bamstreamingmarkduplicates_opts tmpfile=$tmpdir O=" . $markdup_bam_file->path . " index=1 indexfilename=" . $markdup_bai_file->path . " M=" . $markdup_metrics_file->path;
            my $this_cmd       = "use VRPipe::Steps::cram_merge_and_streaming_mark_duplicates; VRPipe::Steps::cram_merge_and_streaming_mark_duplicates->merge_markdup_and_check(q[$preprocess_cmd | $markdup_cmd], file_list_id => $file_list);";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$markdup_bam_file, $markdup_bai_file, $markdup_metrics_file] });
        };
    }
    
    method outputs_definition {
        return {
            markdup_bam_files     => VRPipe::StepIODefinition->create(type => 'bam', max_files => 1, description => 'a BAM file with duplicates marked'),
            markdup_bai_files     => VRPipe::StepIODefinition->create(type => 'bai', max_files => 1, description => 'BAI index file for merged, duplicates marked BAM file'),
            markdup_metrics_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'a text file with metrics from biobambam duplicate marking'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merges CRAM files using samtools merge and marks duplicates using biobambam bamstreamingmarkduplicates";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method merge_markdup_and_check (ClassName|Object $self: Str $cmd_line, Int :$file_list_id! ) {
        my ($out_path) = $cmd_line =~ /O=(\S+)/;
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $out_file = VRPipe::File->get(path => $out_path);
        my @in_files = VRPipe::FileList->get(id => $file_list_id)->files;
        
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        
        my $expected_reads  = 0;
        my $expected_bases  = 0;
        my $expected_paired = 0;
        my $paired          = 0;
        foreach my $in_file (@in_files) {
            my $meta = $in_file->metadata;
            $in_file->disconnect;
            $expected_reads += $meta->{reads} || $in_file->num_records;
            $expected_bases += $meta->{bases} if $meta->{bases};
            $paired = 1 if (exists $meta->{paired});
            $expected_paired ||= $meta->{paired};
        }
        my $actual_reads = $out_file->num_records;
        
        if ($actual_reads == $expected_reads) {
            my %new_meta;
            $new_meta{reads}  = $actual_reads;
            $new_meta{bases}  = $expected_bases if $expected_bases;
            $new_meta{paired} = $expected_paired if $paired;
            $out_file->add_metadata(\%new_meta);
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output CRAM file, yet there were $expected_reads reads in the original CRAM file");
        }
    }
}

1;
