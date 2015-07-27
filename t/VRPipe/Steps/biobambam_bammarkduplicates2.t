#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use POSIX qw(getgroups);

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES SAMTOOLS HTSLIB)],
        required_exe => [qw(bammarkduplicates2 samtools htsfile)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::biobambam_bammarkduplicates2');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('biobambam_bammarkduplicates2', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Marks duplicates for BAM or CRAM files using biobambam bammarkduplicates2'], 'biobambam_bammarkduplicates2 step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'biobambam_bammarkduplicates2',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn',
        method  => 'all',
        source  => file(qw(t data bams.fofn))->absolute->stringify,
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);

ok handle_pipeline(), 'samtools_merge_and_streaming_mark_duplicates pipeline ran ok';

finish;
