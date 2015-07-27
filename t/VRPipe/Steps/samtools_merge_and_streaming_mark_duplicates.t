#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use POSIX qw(getgroups);

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES SAMTOOLS HTSLIB)],
        required_exe => [qw(bamstreamingmarkduplicates samtools htsfile)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::samtools_merge_and_streaming_mark_duplicates');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('samtools_merge_and_streaming_mark_duplicates', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Merges BAM and CRAM files using samtools merge and marks duplicates using biobambam bamstreamingmarkduplicates'], 'samtools_merge_and_streaming_mark_duplicates step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'samtools_merge_and_streaming_mark_duplicates',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => file(qw(t data cram.fofn_with_metadata))->absolute->stringify,
        options => {
            metadata_keys => 'sample',
        }
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);

ok handle_pipeline(), 'samtools_merge_and_streaming_mark_duplicates pipeline ran ok';

$setup = VRPipe::PipelineSetup->create(
    name       => 'samtools_merge_and_streaming_mark_duplicates_on_region',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => file(qw(t data cram.fofn_with_metadata))->absolute->stringify,
        options => {
            metadata_keys => 'sample',
        }
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => { samtools_merge_options => '-R 1:1-10001' }
);

ok handle_pipeline(), 'samtools_merge_and_streaming_mark_duplicates pipeline ran ok on a specified region';

finish;
