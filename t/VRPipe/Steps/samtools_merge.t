#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use POSIX qw(getgroups);

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES SAMTOOLS HTSLIB)],
        required_exe => [qw(samtools htsfile)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::samtools_merge');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('samtools_merge', 'aln_files');
is_deeply [$step->id, $step->description], [1, 'Merges CRAM and BAM files using samtools merge to create a merged BAM file'], 'samtools_merge step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'samtools_merge',
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

ok handle_pipeline(), 'samtools_merge pipeline ran ok';

$setup = VRPipe::PipelineSetup->create(
    name       => 'samtools_merge_on_region',
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
    options     => { samtools_merge_options => '-R 1:1-100' }
);

ok handle_pipeline(), 'samtools_merge pipeline ran ok on a specified region';

finish;
