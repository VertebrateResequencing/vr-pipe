#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 16;
    
    use_ok('VRPipe::Persistent::Schema');
    use_ok('VRPipe::Steps::bam_split_by_sequence');
    
    use TestPipelines;
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_split_by_sequence', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Splits a BAM file into multiple BAM files, one for each sequence the reads were aligned to'], 'bam_split_by_sequence step created and has correct description';

# test using the class methods directly
my $parse_bam = file(qw(t data parser.1kg_lane.bam))->absolute;
my $test_bam = file(qw(t data 2822_6.pe.bam))->absolute;

my $split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $output_dir, pretend => 1);
is scalar keys %$split_bams, 26, 'output number of bams with defaults correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $output_dir, non_chrom => 0, pretend => 1);
is scalar keys %$split_bams, 84, 'output number of bams with non_chrom option unset correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $output_dir, make_unmapped => 1, pretend => 1);
is scalar keys %$split_bams, 27, 'output number of bams with make_unmapped option set correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $output_dir, only => '^20$', pretend => 1);
is scalar keys %$split_bams, 1, 'output number of bams with only chromosome 20 option set correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $output_dir, merge => { '^\d{1,2}$' => 'autosome' }, pretend => 1);
is scalar keys %$split_bams, 5, 'output number of bams with merging autosomal sequence option set correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $output_dir, merge => { '^[^\*]+' => 'mapped' }, non_chrom => 0, make_unmapped => 1, pretend => 1);
is scalar keys %$split_bams, 2, 'output number of bams with merge mapped sequence option set correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $output_dir, ignore => '^\d{1,2}$', pretend => 1);
is scalar keys %$split_bams, 4, 'output number of bams with ignore autosomal sequence option set correct';

my $test_dir = dir($output_dir, 'test_1');
$test_dir->mkpath;
my $ok = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($test_bam, split_dir => $test_dir);
is $ok, 1,  'split_bam_by_sequence() default test ran okay';

$test_dir = dir($output_dir, 'test_2');
$test_dir->mkpath;
$ok = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($test_bam, split_dir => $test_dir, make_unmapped => 1);
is $ok, 1,  'split_bam_by_sequence() make_unmapped ran okay';

$test_dir = dir($output_dir, 'test_3');
$test_dir->mkpath;
$ok = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($test_bam, split_dir => $test_dir, ignore => 'chr2');
is $ok, 1,  'split_bam_by_sequence() ignore test ran okay';

$test_dir = dir($output_dir, 'test_4');
$test_dir->mkpath;
$ok = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($test_bam, split_dir => $test_dir, only => 'chr2');
is $ok, 1,  'split_bam_by_sequence() only test ran okay';

$test_dir = dir($output_dir, 'test_5');
$test_dir->mkpath;
$ok = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($test_bam, split_dir => $test_dir, merge => { '^[^\*]+' => 'mapped' }, make_unmapped => 1);
is $ok, 1,  'split_bam_by_sequence() merge test ran okay';

# test as part of a pipeline
my $setup = VRPipe::PipelineSetup->get(name => 'split_setup',
                                       datasource => VRPipe::DataSource->get(type => 'fofn', method => 'all', source => file(qw(t data improvement_datasource.fofn))->absolute),
                                       output_root => $output_dir,
                                       pipeline => $pipeline,
                                       options => { split_bam_merge => q[{ '^[^\\*]+' => 'mapped' }], split_bam_make_unmapped => 1 });

my (@output_files);
my $element_id=0;
foreach my $bam (qw(2822_7.pe.bam 2822_6.pe.bam 2822_6.se.bam 2823_4.pe.bam 8324_8.pe.bam)) {
    $element_id++;
    foreach my $prefix (qw(mapped unmapped)) {
        push(@output_files, file($output_dir, output_subdirs($element_id), '1_bam_split_by_sequence', "$prefix.$bam"));
    }
}
ok handle_pipeline(@output_files), 'single step bam_split_by_sequence pipeline ran and created all expected output files';

finish;
