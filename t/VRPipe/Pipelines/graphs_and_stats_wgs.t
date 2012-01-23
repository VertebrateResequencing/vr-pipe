#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => 'VRPIPE_TEST_PIPELINES');
    use TestPipelines;
}

#./Build test --test_files t/VRPipe/Pipelines/graphs_and_stats_wgs.t --verbose

my $stats_output_dir = get_output_dir('graphs_and_stats_wgs_pipeline');

ok my $stats_pipeline = VRPipe::Pipeline->get(name => 'graphs_and_stats_wgs'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($stats_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(qc_stats_bamcheck_wgs
							 qc_stats_bamcheck_rmdup_wgs
							 qc_plots_wgs);

is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data human_g1k_v37.chr20.fa));
my $ref_dir = dir($stats_output_dir, 'ref');
$stats_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'human_g1k_v37.chr20.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $stats_fa_source = file(qw(t data human_g1k_v37.chr20.fa.stats));
my $ref_fa_stats = file($ref_dir, 'human_g1k_v37.chr20.fa.stats')->stringify;
copy($stats_fa_source, $ref_fa_stats);

my $fa_index_source = file(qw(t data human_g1k_v37.chr20.fa.fai));
my $fa_index = file($ref_dir, 'human_g1k_v37.chr20.fa.fai')->stringify;
copy($fa_index_source, $fa_index);

my $stats_pipelinesetup = VRPipe::PipelineSetup->get(name => 'graphs_and_stats',
                                                       datasource => VRPipe::DataSource->get(type => 'fofn',
                                                       method => 'all',
                                                       source => file(qw(t data hs_chr20.qc.bam.fofn))->absolute),
                                                       output_root => $stats_output_dir,
                                                       pipeline => $stats_pipeline,
                                                       options => {reference_fasta => $ref_fa, reference_fasta_stats =>$ref_fa_stats});

my (@output_files,@final_files);
my $element_id=0;
my $sample = "NA20526";

foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
	$element_id++;
    push(@output_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "${in}.bam.bamcheck"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-coverage.gp"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-gc-content.png"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-insert-size.gp"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-quals3.png"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-quals-hm.png"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-acgt-cycles.gp"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-coverage.png"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-gc-depth.gp"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-insert-size.png"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-quals.gp"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-quals.png"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-acgt-cycles.png"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-gc-content.gp"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-gc-depth.png"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-quals2.png"));
    push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_stats_bamcheck_rmdup_wgs', "$sample-quals-hm.gp"));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;

finish;