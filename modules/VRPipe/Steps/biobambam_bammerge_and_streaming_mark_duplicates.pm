
=head1 NAME

VRPipe::Steps::biobambam_bammerge_and_streaming_mark_duplicates - a step

=head1 DESCRIPTION

Runs biobambam's bammerge and bamstreamingmarkduplicates to mark duplicates.

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

class VRPipe::Steps::biobambam_bammerge_and_streaming_mark_duplicates with VRPipe::StepRole {
    method options_definition {
        return {
            bammerge_exe                       => VRPipe::StepOption->create(description => 'path to bammerge executable',                                                               optional => 1, default_value => 'bammerge'),
            bammerge_options                   => VRPipe::StepOption->create(description => 'options for bammerge',                                                                      optional => 1, default_value => 'level=0'),
            bamstreamingmarkduplicates_exe     => VRPipe::StepOption->create(description => 'path to bamstreamingmarkduplicates executable',                                             optional => 1, default_value => 'bamstreamingmarkduplicates'),
            bamstreamingmarkduplicates_options => VRPipe::StepOption->create(description => 'bamstreamingmarkduplicates options (excluding arguments that set input/output file names)', optional => 1, default_value => 'resetdupflag=1 outputthreads=4'),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'aln',                                          # cram or bam
                max_files   => -1,
                description => '1 or more coordinate sorted BAM or CRAM files',
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
            my $bammerge                        = $options->{bammerge_exe};
            my $bammerge_opts                   = $options->{bammerge_options};
            my $bamstreamingmarkduplicates      = $options->{bamstreamingmarkduplicates_exe};
            my $bamstreamingmarkduplicates_opts = $options->{bamstreamingmarkduplicates_options};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bammerge',
                    version => VRPipe::StepCmdSummary->determine_version($bammerge . ' --version', '^This is biobambam\d? version (.+)\.$'),
                    summary => "bammerge $bammerge_opts \$input_file(s)"
                )
            );
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bamstreamingmarkduplicates',
                    version => VRPipe::StepCmdSummary->determine_version($bamstreamingmarkduplicates . ' --version', '^This is biobambam\d? version (.+)\.$'),
                    summary => "bamstreamingmarkduplicates $bamstreamingmarkduplicates_opts O=\$markdup_file index=1 indexfilename=\$markdup_index M=\$metrics_file"
                )
            );
            
            my $cpus             = 1;
            my ($input_threads)  = $bammerge_opts =~ m/inputthreads=(\d+)/;
            my ($output_threads) = $bammerge_opts =~ m/outputthreads=(\d+)/;
            $cpus = $input_threads  if ($input_threads  && $input_threads > $cpus);
            $cpus = $output_threads if ($output_threads && $output_threads > $cpus);
            ($input_threads)  = $bamstreamingmarkduplicates_opts =~ m/inputthreads=(\d+)/;
            ($output_threads) = $bamstreamingmarkduplicates_opts =~ m/outputthreads=(\d+)/;
            $cpus = $input_threads  if ($input_threads  && $input_threads > $cpus);
            $cpus = $output_threads if ($output_threads && $output_threads > $cpus);
            my $req = $self->new_requirements(memory => 8000, time => 1, cpus => $cpus);
            
            my $inputs          = $self->inputs->{bam_files};
            my $merged_metadata = $self->common_metadata($inputs);
            my @input_paths     = map { "I=" . $_->path } @$inputs;
            my $file_list       = VRPipe::FileList->create(files => $inputs)->id;
            my $basename        = $self->step_state->dataelement->id;
            my $markdup_file    = $self->output_file(
                output_key => 'markdup_files',
                basename   => $basename . '.bam',
                type       => 'bam',
                metadata   => $merged_metadata
            );
            my $markdup_index_file = $self->output_file(
                output_key => 'markdup_index_files',
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
            my $markdup_path = $markdup_file->path;
            my $merge_cmd    = qq[$bammerge $bammerge_opts @input_paths];
            my $markdup_cmd  = "$bamstreamingmarkduplicates $bamstreamingmarkduplicates_opts O=$markdup_path index=1 indexfilename=" . $markdup_index_file->path . " M=" . $markdup_metrics_file->path;
            my $this_cmd     = "use VRPipe::Steps::biobambam_bammerge_and_streaming_mark_duplicates; VRPipe::Steps::biobambam_bammerge_and_streaming_mark_duplicates->merge_markdup_and_check(q[$merge_cmd | $markdup_cmd], input_file_list => $file_list, out_file => q[$markdup_path]);";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$markdup_file, $markdup_index_file, $markdup_metrics_file] });
        };
    }
    
    method outputs_definition {
        return {
            markdup_files         => VRPipe::StepIODefinition->create(type => 'aln', max_files => 1, description => 'a merged BAM file with duplicates marked'),
            markdup_index_files   => VRPipe::StepIODefinition->create(type => 'bin', max_files => 1, description => ' index file for merged, markdup-ed BAM file'),
            markdup_metrics_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'a text file with metrics from biobambam duplicate marking'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merges BAM and CRAM files using biobambam bammerge and marks duplicates using biobambam bamstreamingmarkduplicates";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method merge_markdup_and_check (ClassName|Object $self: Str $cmd_line, Int :$input_file_list!, Str|File :$out_file! ) {
        unless (ref($out_file) && ref($out_file) eq 'VRPipe::File') {
            $out_file = VRPipe::File->create(path => file($out_file));
        }
        
        my $out_path  = $out_file->path;
        my $out_index = VRPipe::File->get(path => $out_file->type eq 'bam' ? qq[$out_path.bai] : qq[$out_path.crai]);
        my @in_files  = VRPipe::FileList->get(id => $input_file_list)->files;
        
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        $out_index->update_stats_from_disc(retries => 3);
        
        my $actual_reads    = $out_file->num_records;
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
