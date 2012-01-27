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

#./Build test --test_files t/VRPipe/Pipelines/graphs_and_stats_exome.t --verbose

my $stats_output_dir = get_output_dir('graphs_and_stats_exome_pipeline');

ok my $stats_pipeline = VRPipe::Pipeline->get(name => 'graphs_and_stats_exome'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($stats_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(qc_stats_exome
							 qc_plots_exome);

is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $intervals_dump = file(qw(t data exome.qc.chr20.dump))->absolute->stringify;

my $stats_pipelinesetup = VRPipe::PipelineSetup->get(name => 'graphs_and_stats_exome',
                                                       datasource => VRPipe::DataSource->get(type => 'fofn',
                                                                                             method => 'all',
                                                                                             source => file(qw(t data hs_chr20.qc.bam.fofn))->absolute),
                                                       output_root => $stats_output_dir,
                                                       pipeline => $stats_pipeline,
                                                       options => { load_intervals_dump_file => $intervals_dump });

my (@output_files,@final_files);
my $element_id=0;
my $sample = "NA20526";

foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
	$element_id++;
    push(@output_files, file($stats_output_dir, output_subdirs($element_id), '1_qc_stats_exome', "${in}.bam.stats.dump"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.bait_gc_vs_cvg.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.bait_gc_vs_cvg.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.bait_gc_vs_cvg.scaled.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.bait_gc_vs_cvg.scaled.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.coverage_per_base.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.coverage_per_base.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.cumulative_coverage.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.cumulative_coverage.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.gc_mapped.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.gc_mapped.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.gc_unmapped.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.gc_unmapped.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.insert_size.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.insert_size.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.mean_coverage.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.mean_coverage.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.normalised_coverage.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.normalised_coverage.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.quality_scores_1.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.quality_scores_1.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.quality_scores_2.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.quality_scores_2.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.target_gc_vs_cvg.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.target_gc_vs_cvg.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.target_gc_vs_cvg.scaled.png"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "$sample.target_gc_vs_cvg.scaled.png.R"));
	push(@final_files, file($stats_output_dir, output_subdirs($element_id), '2_qc_plots_exome', "Rplots.pdf"));	
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;

finish;
