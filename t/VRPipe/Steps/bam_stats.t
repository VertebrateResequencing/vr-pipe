#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 12;
    use VRPipeTest (required_env => [qw(SAMTOOLS)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::bam_stats');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_stats', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Calculates various statistics about bam files, producing .bas files'], 'bam_stats step created and has correct description';

# test using the class methods directly
my $test_bam = file(qw(t data bas.bam))->absolute;
is_deeply { VRPipe::Steps::bam_stats->bam_statistics($test_bam) },
  {
    SRR00001 => {
        total_bases                  => 115000,
        mapped_bases                 => 58583,
        total_reads                  => 2000,
        mapped_reads                 => 1084,
        mapped_reads_paired_in_seq   => 1084,
        mapped_reads_properly_paired => 1070,
        percent_mismatch             => '2.05',
        avg_qual                     => '23.32',
        avg_isize                    => 286,
        sd_isize                     => '74.10',
        median_isize                 => 275,
        mad                          => 48,
        duplicate_reads              => 2,
        duplicate_bases              => 122
    }
  },
  'bam_statistics test';

my $given_bas = file($output_dir, 'test.bas');
my $ok = VRPipe::Steps::bam_stats->bas($test_bam, $given_bas, release_date => 20100208);
is $ok, 1, 'bas ran ok';
my $expected_bas = file(qw(t data example.bas));
ok open(my $ebfh, $expected_bas), 'opened expected .bas';
my @expected = <$ebfh>;
close($ebfh);
ok open(my $tbfh, $given_bas), 'opened result .bas';
my @given = <$tbfh>;
close($tbfh);
is_deeply \@given, \@expected, 'bas output was as expected';

# test making a bas file from a bam with RG in PU instead of ID
$given_bas = file($output_dir, 'test2.bas');
$ok = VRPipe::Steps::bam_stats->bas(file(qw(t data rg_pu.bam))->absolute, $given_bas, release_date => 20110521, rg_from_pu => 1);
is $ok, 1, 'bas() with rg_from_pu ran ok';
$expected_bas = file(qw(t data example3.bas));
ok open($ebfh, $expected_bas), 'opened expected .bas';
my @expected2 = <$ebfh>;
close($ebfh);
ok open($tbfh, $given_bas), 'opened result .bas';
@given = <$tbfh>;
close($tbfh);
is_deeply \@given, \@expected2, 'bas output was as expected when converting RG PU to ID';

SKIP: {
    my $num_tests = 1;
    skip "longer-running pipeline test disabled without VRPIPE_TEST_PIPELINES", $num_tests unless $ENV{VRPIPE_TEST_PIPELINES};
    
    # test as part of a pipeline
    my $setup = VRPipe::PipelineSetup->create(
        name        => 'bsp_setup',
        datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data datasource.bam_fofn))->absolute),
        output_root => $output_dir,
        pipeline    => $pipeline,
        options     => {}
    );
    
    ok handle_pipeline(), 'single-step pipeline ran ok';
}

finish;
