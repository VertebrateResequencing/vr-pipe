#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 8;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(gunzip)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::fastq_decompress');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('fastq_decompress', 'compressed_fastq_files');
is_deeply [$step->id, $step->description], [1, 'Decompresses fastq files'], 'fastq_decompress step created and has correct description';

# test using the class methods directly
my $test_fq = VRPipe::File->create(path => file(qw(t data 2822_6.fastq.gz))->absolute)->path;
my $out_fq = VRPipe::File->create(path => file($output_dir, '2822_6.fastq')->absolute)->path;

my $cmd = "gunzip -c $test_fq > $out_fq";

my $ok = VRPipe::Steps::fastq_decompress->decompress_and_check($cmd);
is $ok, 1, 'decompress_and_check ran ok';

my $in_file = VRPipe::File->create(path => $test_fq);
ok $in_file, 'got test fq file object';
is $in_file->num_records, 50, 'expected number of records in test fq';

my $out_file = VRPipe::File->create(path => $out_fq);
ok $out_file, 'got output fq file object';
is $out_file->num_records, 50, 'expected number of records in output fq';

# test as part of a pipeline
my $setup = VRPipe::PipelineSetup->create(
    name        => 'fq_setup',
    datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data decompress_datasource.fofn))->absolute),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);

ok handle_pipeline(), 'single-step pipeline ran ok';

finish;
