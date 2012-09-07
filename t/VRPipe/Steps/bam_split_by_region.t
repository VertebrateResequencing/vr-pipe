#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES SAMTOOLS)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::bam_split_by_region');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_split_by_region', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Splits a BAM file into multiple BAM files, one for each genomic region specified.'], 'bam_split_by_region step created and has correct description';

# test using the class methods directly
my $test_bam  = VRPipe::File->create(path => file(qw(t data NA19334.bam))->absolute)->path->stringify;
my $ref_index = VRPipe::File->create(path => file(qw(t data human_g1k_v37.chr11.chr20.fa.gz.fai))->absolute)->path->stringify;
my $split_dir = $output_dir->stringify;

my $ok = VRPipe::Steps::bam_split_by_region->split_bam_by_region($test_bam, split_dir => $split_dir, include_mate => 1, regions => ['11_1-100000', '11_98000-200000']);
is $ok, 1, 'split_bam_by_region() default test ran okay';

# test as part of a pipeline
VRPipe::PipelineSetup->create(
    name        => 'split_by_region_step',
    datasource  => VRPipe::DataSource->create(type => 'fofn_with_metadata', method => 'all', source => file(qw(t data calling_datasource.fofn))->absolute),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_index                   => $ref_index,
        split_bam_by_region_chunk_size    => 10000000,
        split_bam_by_region_chunk_overlap => 5000,
        split_bam_by_region_chrom_list    => '11 20',
        include_mate                      => 1
    }
);

my $regions = [qw(11_1-10000000 11_9995001-19995000 11_19990001-29990000 11_29985001-39985000 11_39980001-49980000 11_49975001-59975000 11_59970001-69970000 11_69965001-79965000 11_79960001-89960000 11_89955001-99955000 11_99950001-109950000 11_109945001-119945000 11_119940001-129940000 11_129935001-135006516 20_1-10000000 20_19990001-29990000 20_29985001-39985000 20_39980001-49980000 20_49975001-59975000 20_59970001-63025520 20_9995001-19995000)];

my (@output_files);
my $element_id = 0;
foreach my $sample (qw(NA19334 NA19381)) {
    $element_id++;
    foreach my $region (@$regions) {
        push(@output_files, file(output_subdirs($element_id), '1_bam_split_by_region', "$region.$sample.bam"));
    }
}
ok handle_pipeline(@output_files), 'single step bam_split_by_region pipeline ran and created all expected output files';

finish;
