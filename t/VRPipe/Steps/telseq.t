#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use POSIX qw(getgroups);

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(telseq)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::telseq');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('telseq', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Runs telseq on BAM files to estimate telomere length'], 'telseq step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'telseq test',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'all',
        source  => file(qw(t data calling_datasource.fofn))->absolute->stringify,
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);

ok handle_pipeline(), 'telseq pipeline ran ok';

finish;
