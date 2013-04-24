#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES CONVEX_R_LIB)],
        max_retries  => 1,
        required_exe => [qw(R)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('convex_plot_generation_pipeline');

ok my $pipeline1 = VRPipe::Pipeline->create(name => 'convex_plot_generation'), 'able to create the convex_plot_generation pipeline';

my @s_names;
foreach my $stepmember ($pipeline1->steps) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(convex_plots);
is_deeply \@s_names, \@expected_step_names, 'the rd pipeline has the correct steps';

my $convex_r_libs = $ENV{CONVEX_R_LIB};

my $fofn = file(qw(t data cnv cnv.calls.fofn))->absolute;

my $pipelinesetup1 = VRPipe::PipelineSetup->create(
    name        => 'convex_plot_generation_pipeline',
    datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'group_all', source => $fofn),
    output_root => $output_dir,
    pipeline    => $pipeline1,
    options     => {
        cleanup => 0,
        r_libs  => "$convex_r_libs",
    }
);

my (@output_files, @final_files);
my @output_subdirs = output_subdirs(1);
push(@output_files, file(@output_subdirs, '1_convex_plots', 'CnvFofn.txt'));
push(@output_files, file(@output_subdirs, '1_convex_plots', 'CNVstats_CallsperSample.png'));
push(@output_files, file(@output_subdirs, '1_convex_plots', 'CNVstats_DelDupRatio.png'));

ok handle_pipeline(@output_files, @final_files), 'rd pipeline ran and created all expected output files';

done_testing;
