#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Cwd;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('pluritest_gene_expression_analysis_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'pluritest_gene_expression_analysis'), 'able to get the pluritest_gene_expression_analysis pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(genome_studio_expression_reformat plot_pluritest_gene_expression);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

ok my $ds = VRPipe::DataSource->create(
    type   => 'fofn',
    method => 'all',
    source => file(qw(t data hipsci_gen_ex.fofn))->absolute
), 'could create a fofn datasource';

#issues with R on precise-dev64 
#  -> path to working R is /software/bin/R-2.15.2, but this does not exist on uk10k
#  -> no lumi package in this version ?
my $r_bin_path = '/software/bin/R-2.15';
my $pluritest_script = '/software/vertres/scripts/pluriTest_commandLine_vrpipe.r';

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my pluritest_gene_expression_analysis_pipeline',
    datasource => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        'reformat_annotation'    => file(qw(t data hipsci_gene_expression_annotation.txt))->absolute->stringify,
        'reformat_mapping'       => file(qw(t data hipsci_gene_expression_mapping.txt))->absolute->stringify,
        'pluritest_script'       => $pluritest_script,
        'pluritest_data'         => file(qw(t data pluritest.RData))->absolute->stringify,
        'r_bin_path'             => $r_bin_path,
        cleanup                  => 0
    }
);

my (@output_files, @final_files);
my @output_dirs = output_subdirs(1);
push(@output_files, file(@output_dirs, '1_genome_studio_expression_reformat', "hipsci_gene_expression_profile.reformat.txt"));
foreach my $kind (qw(01 02a 02 03c 03)) {
	push(@final_files, file(@output_dirs, '2_plot_pluritest_gene_expression', 'pluritest_image' . $kind . '.png'));
}
push(@final_files,  file(@output_dirs, '2_plot_pluritest_gene_expression', 'pluritest.csv'));

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
