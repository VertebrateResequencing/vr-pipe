
=head1 NAME

VRPipe::Steps::samtools_merge_and_streaming_mark_duplicates - a step

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

class VRPipe::Steps::samtools_merge_and_streaming_mark_duplicates extends VRPipe::Steps::samtools_merge {
    around options_definition {
        return {
            %{ $self->$orig },
            samtools_merge_options          => VRPipe::StepOption->create(description => 'options for samtools merge',                                                                optional => 1, default_value => '-u'),
            bamstreamingmarkduplicates_exe  => VRPipe::StepOption->create(description => 'path to bamstreamingmarkduplicates executable',                                             optional => 1, default_value => 'bamstreamingmarkduplicates'),
            bamstreamingmarkduplicates_opts => VRPipe::StepOption->create(description => 'bamstreamingmarkduplicates options (excluding arguments that set input/output file names)', optional => 1, default_value => 'resetdupflag=1'),
        };
    }
    
    method body_sub {
        return sub {
            my $self                            = shift;
            my $options                         = $self->options;
            my $samtools                        = $options->{samtools_exe};
            my $samtools_merge_opts             = $options->{samtools_merge_options};
            my $bamstreamingmarkduplicates_exe  = $options->{bamstreamingmarkduplicates_exe};
            my $bamstreamingmarkduplicates_opts = $options->{bamstreamingmarkduplicates_opts};
            my $check_read_sum                  = $samtools_merge_opts =~ m/-R/ ? 0 : 1;
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => "samtools merge -u $samtools_merge_opts - \$input_file(s)"
                )
            );
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bamstreamingmarkduplicates',
                    version => VRPipe::StepCmdSummary->determine_version($bamstreamingmarkduplicates_exe . ' --version', '^This is biobambam\d? version (.+)\.$'),
                    summary => "bamstreamingmarkduplicates $bamstreamingmarkduplicates_opts O=\$markdup_bam index=1 indexfilename=\$markdup_bam.bai M=\$metrics_file"
                )
            );
            
            my $req = $self->new_requirements(memory => 3000, time => 1);
            
            my $inputs          = $self->inputs->{aln_files};
            my $merged_metadata = $self->common_metadata($inputs);
            $merged_metadata = { %$merged_metadata, $self->element_meta };
            my @input_paths = map { $_->path } @$inputs;
            my $file_list = VRPipe::FileList->create(files => $inputs)->id;
            my $basename = $self->step_state->dataelement->id;
            if (defined $$merged_metadata{chrom} && defined $$merged_metadata{from} && defined $$merged_metadata{to}) {
                my ($chrom, $from, $to) = ($$merged_metadata{chrom}, $$merged_metadata{from}, $$merged_metadata{to});
                $samtools_merge_opts .= " -R $chrom:$from-$to";
                $basename       = "${chrom}_${from}-${to}.$basename";
                $check_read_sum = 0;
            }
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
            my $markdup_path = $markdup_bam_file->path;
            my $merge_cmd    = qq[$samtools merge $samtools_merge_opts - @input_paths];
            my $markdup_cmd  = "$bamstreamingmarkduplicates_exe $bamstreamingmarkduplicates_opts O=$markdup_path index=1 indexfilename=" . $markdup_bai_file->path . " M=" . $markdup_metrics_file->path;
            my $this_cmd     = "use VRPipe::Steps::samtools_merge_and_streaming_mark_duplicates; VRPipe::Steps::samtools_merge_and_streaming_mark_duplicates->merge_and_check(q[$merge_cmd | $markdup_cmd], input_file_list => $file_list, out_file => q[$markdup_path], check_read_sum => $check_read_sum);";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$markdup_bam_file, $markdup_bai_file, $markdup_metrics_file] });
        };
    }
    
    method outputs_definition {
        return {
            markdup_bam_files     => VRPipe::StepIODefinition->create(type => 'bam', max_files => 1, description => 'a BAM file with duplicates marked'),
            markdup_bai_files     => VRPipe::StepIODefinition->create(type => 'bai', max_files => 1, description => 'BAI index file for merged, markdup-ed BAM file'),
            markdup_metrics_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'a text file with metrics from biobambam duplicate marking'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merges BAM and CRAM files using samtools merge and marks duplicates using biobambam bamstreamingmarkduplicates";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }

}

1;
