
=head1 NAME

VRPipe::Steps::biobambam_mark_duplicates_two - a step

=head1 DESCRIPTION

Runs the biobambam's mark duplicates 2 to mark duplicates.

=head1 AUTHOR

Martin Pollard <mp15@sanger.ac.uk>.

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

class VRPipe::Steps::biobambam_mark_duplicates_two with VRPipe::StepRole {

    method options_definition {
        return {
            bammarkduplicates2_exe     => VRPipe::StepOption->create(description => 'path to bammarkduplicates2 executable', optional => 1, default_value => 'bammarkduplicates2'),
            bammarkduplicates2_opts    => VRPipe::StepOption->create(description => 'bammarkduplicates2 options (excluding arguments that set input/output file names and threads)', optional => 1, default_value => ''),
            bammarkduplicates2_threads => VRPipe::StepOption->create(description => 'number of threads to use', optional => 1, default_value => '1'),
            tmp_dir                    => VRPipe::StepOption->create(description => 'location for tmp directories; defaults to working directory', optional => 1),
        };
    }

    method inputs_definition {
        return {
            cram_files => VRPipe::StepIODefinition->create(
                type        => 'cram',
                max_files   => -1,
                description => '1 or more coordinate sorted CRAM files',
                metadata    => {
                    lane             => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    bases            => 'total number of base pairs',
                    reads            => 'total number of reads (sequences)',
                    forward_reads    => 'number of forward reads',
                    reverse_reads    => 'number of reverse reads',
                    paired           => '0=single ended reads only; 1=paired end reads present',
                    mean_insert_size => 'mean insert size (0 if unpaired)',
                    avg_read_length  => 'the average length of reads',
                    library          => 'library name',
                    sample           => 'sample name',
                    center_name      => 'center name',
                    platform         => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study            => 'name of the study, put in the DS field of the RG header line',
                    optional         => ['lane', 'library', 'sample', 'center_name', 'platform', 'study', 'mean_insert_size', 'forward_reads', 'reverse_reads', 'paired', 'avg_read_length']
                }
            ),
        };
    }

    method body_sub {
        return sub {
            my $self            = shift;
            my $options         = $self->options;
            my $bammarkduplicates2_exe     = $options->{bammarkduplicates2_exe};
            my $bammarkduplicates2_opts    = $options->{bammarkduplicates2_opts};
            my $bammarkduplicates2_threads = $options->{bammarkduplicates2_threads};

            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bammarkduplicates2',
                    version => VRPipe::StepCmdSummary->determine_version($bammarkduplicates2_exe . ' --version', '^This is biobambam version (.+)\.$'),
                    summary => "bammarkduplicates2 $bammarkduplicates2_opts I:\$bam_file O:\$output_file(s)"
                )
            );

            my $req = $self->new_requirements(memory => 3000, time => 1, cpus => $bammarkduplicates2_threads);

            my $input_files = join ' ', (map {$_ = "I=".$_->path;} (@{ $self->inputs->{cram_files} }));
            my $merged_metadata = $self->common_metadata($self->inputs->{cram_files});
            my $basename = $self->stepstate->dataelement->id . '.cram';
            my $markdup_bam_file = $self->output_file(
                output_key => 'markdup_cram_files',
                basename   => 'mark_merge',
                type       => 'cram',
                metadata   => $merged_metadata
            );
            my $tmpdir = $options->{tmp_dir} || $markdup_bam_file->dir;
            my $threads = "";
            $threads = " markthreads=$bammarkduplicates2_threads" if $bammarkduplicates2_threads ne '1';
            my $this_cmd = "$bammarkduplicates2_exe ${bammarkduplicates2_opts}${threads} tmpfile=$tmpdir $input_files O=" . $markdup_bam_file->path . " $bammarkduplicates2_opts";

            $self->dispatch_wrapped_cmd('VRPipe::Steps::biobambam_mark_duplicates_two', 'markdup_and_check', [$this_cmd, $req, { output_files => [$markdup_bam_file] }]);
        };
    }

    method outputs_definition {
        return {
            markdup_cram_files => VRPipe::StepIODefinition->create(type => 'cram', max_files => 1, description => 'a cram file with duplicates marked'),
        };
    }

    method post_process_sub {
        return sub { return 1; };
    }

    method description {
        return "Marks duplicates in CRAM files using biobambam";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }

    method markdup_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($out_path) = $cmd_line =~ /O=(\S+)/;
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        my @inputs = $cmd_line =~ /I=(\S+)/g;
        (scalar @inputs)  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        my $expected_reads = 0;
        foreach my $in_path (@inputs) {
            my $in_file  = VRPipe::File->get(path => $in_path);
            $in_file->disconnect;
            $expected_reads += $in_file->metadata->{reads} || $in_file->num_records;
        }
        my $out_file = VRPipe::File->get(path => $out_path);

        system($cmd_line) && $self->throw("failed to run [$cmd_line]");

        $out_file->update_stats_from_disc(retries => 3);
        my $actual_reads = $out_file->num_records;

        if ($actual_reads == $expected_reads) {
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output CRAM file, yet there were $expected_reads reads in the original CRAM file");
        }
    }
}

1;
