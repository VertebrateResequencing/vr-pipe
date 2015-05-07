#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use POSIX qw(getgroups);

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES SAMTOOLS)],
        required_exe => [qw(bamstreamingmarkduplicates samtools)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::cram_merge_and_streaming_mark_duplicates');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('cram_merge_and_streaming_mark_duplicates', 'cram_files');
is_deeply [$step->id, $step->description], [1, 'Merges CRAM files using samtools merge and marks duplicates using biobambam bamstreamingmarkduplicates'], 'cram_merge_and_streaming_mark_duplicates step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'cram_merge_and_streaming_mark_duplicates',
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

ok handle_pipeline(), 'cram_merge_and_streaming_mark_duplicates pipeline ran ok';

finish;
