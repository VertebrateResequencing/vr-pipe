#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 19;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES SAMTOOLS)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::bam_split_by_sequence');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_split_by_sequence', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Splits a BAM file into multiple BAM files, one for each sequence the reads were aligned to'], 'bam_split_by_sequence step created and has correct description';

# test using the class methods directly
my $parse_bam = VRPipe::File->create(path => file(qw(t data parser.1kg_lane.bam))->absolute)->path->stringify;
my $test_bam = VRPipe::File->create(path => file(qw(t data 2822_6.pe.bam))->absolute)->path->stringify;
my $split_dir = $output_dir->stringify;

my $split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $split_dir, pretend => 1);
is scalar keys %$split_bams, 26, 'output number of bams with defaults correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $split_dir, non_chrom => 0, pretend => 1);
is scalar keys %$split_bams, 84, 'output number of bams with non_chrom option unset correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $split_dir, make_unmapped => 1, pretend => 1);
is scalar keys %$split_bams, 27, 'output number of bams with make_unmapped option set correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $split_dir, only => '^20$', pretend => 1);
is scalar keys %$split_bams, 1, 'output number of bams with only chromosome 20 option set correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $split_dir, merge => { '^\d{1,2}$' => 'autosome' }, pretend => 1);
is scalar keys %$split_bams, 5, 'output number of bams with merging autosomal sequence option set correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $split_dir, merge => { '^[^\*]+' => 'mapped' }, non_chrom => 0, make_unmapped => 1, pretend => 1);
is scalar keys %$split_bams, 2, 'output number of bams with merge mapped sequence option set correct';

$split_bams = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($parse_bam, split_dir => $split_dir, ignore => '^\d{1,2}$', pretend => 1);
is scalar keys %$split_bams, 4, 'output number of bams with ignore autosomal sequence option set correct';

my $test_dir = dir($output_dir, 'test_1');
$test_dir->mkpath;
my $ok = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($test_bam, split_dir => $test_dir->stringify);
is $ok, 1,  'split_bam_by_sequence() default test ran okay';

$test_dir = dir($output_dir, 'test_2');
$test_dir->mkpath;
$ok = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($test_bam, split_dir => $test_dir->stringify, make_unmapped => 1);
is $ok, 1,  'split_bam_by_sequence() make_unmapped ran okay';

$test_dir = dir($output_dir, 'test_3');
$test_dir->mkpath;
$ok = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($test_bam, split_dir => $test_dir->stringify, ignore => 'chr2');
is $ok, 1,  'split_bam_by_sequence() ignore test ran okay';

$test_dir = dir($output_dir, 'test_4');
$test_dir->mkpath;
$ok = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($test_bam, split_dir => $test_dir->stringify, only => 'chr2');
is $ok, 1,  'split_bam_by_sequence() only test ran okay';

$test_dir = dir($output_dir, 'test_5');
$test_dir->mkpath;
$ok = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($test_bam, split_dir => $test_dir->stringify, merge => { '^[^\*]+' => 'mapped' }, make_unmapped => 1);
is $ok, 1,  'split_bam_by_sequence() merge test ran okay';

# test as part of a pipeline
VRPipe::PipelineSetup->create(name => 'split_setup',
                           datasource => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data improvement_datasource.fofn))->absolute),
                           output_root => $output_dir,
                           pipeline => $pipeline,
                           options => { split_bam_merge => q[{ '^fake_chr[12]$' => 'autosome', '^[^\\*]+' => 'mapped' }], split_bam_make_unmapped => 1 });

my (@output_files);
my $element_id=0;
foreach my $bam (qw(2822_7.pe.bam 2822_6.pe.bam 2822_6.se.bam 2823_4.pe.bam 8324_8.pe.bam)) {
   $element_id++;
   foreach my $prefix (qw(autosome mapped unmapped)) {
       push(@output_files, file(output_subdirs($element_id), '1_bam_split_by_sequence', "$prefix.$bam"));
   }
}
ok handle_pipeline(@output_files), 'single step bam_split_by_sequence pipeline ran and created all expected output files';


VRPipe::DataElementState->create(pipelinesetup => 1, dataelement => 5)->start_from_scratch();

my @file_exists = map { -s $_ ? 1 : 0 } @output_files;
is_deeply \@file_exists, [1,1,1,1,1,1,1,1,1,1,1,1,0,0,0], 'correct files were removed after start from scratch';

ok handle_pipeline(@output_files), 'all file recreated after an element was started from scratch';

# run a longer test
VRPipe::PipelineSetup->create(name => 'split_setup2',
                           datasource => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data split_bam_datasource.fofn))->absolute),
                           output_root => $output_dir,
                           pipeline => $pipeline,
                           options => { split_bam_merge => q[{'^1$' => 'chrom1', '^2$' => 'chrom2', '^3$'  => 'chrom3', '^4$'  => 'chrom4', '^5$'  => 'chrom5', '^6$'  => 'chrom6', '^7$'  => 'chrom7', '^8$'  => 'chrom8', '^9$'  => 'chrom9', '^10$' => 'chrom10', '^11$' => 'chrom11', '^12$' => 'chrom12', '^13$' => 'chrom13', '^14$' => 'chrom14', '^15$' => 'chrom15', '^16$' => 'chrom16', '^17$' => 'chrom17', '^18$' => 'chrom18', '^19$' => 'chrom19', '^20$' => 'chrom20', '^21$' => 'chrom21', '^22$' => 'chrom22', '^X$'  => 'chromX', '^Y$'  => 'chromY', '^MT$'  => 'chromMT', '^[^\\*]+' => 'mapped' }], 
                                        split_bam_make_unmapped => 1,
                                        split_bam_non_chrom => 0, });

@output_files = ();
$element_id++;
foreach my $chrom (1..22,qw(X Y MT)) {
    push(@output_files, file(output_subdirs($element_id, 2), '1_bam_split_by_sequence', "chrom$chrom.g1k_small.bam"));
}
foreach my $chrom (qw(mapped unmapped)) {
    push(@output_files, file(output_subdirs($element_id, 2), '1_bam_split_by_sequence', "$chrom.g1k_small.bam"));
}
ok handle_pipeline(@output_files), 'single step bam_split_by_sequence pipeline ran and created all expected output files for longer g1k test';

my @bam_records = get_bam_records($output_files[6]);
is_deeply \@bam_records, [], 'header only bam was correctly created when no sequence was present in the bam file';

finish;
