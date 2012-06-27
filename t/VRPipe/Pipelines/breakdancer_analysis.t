#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],
	);
    use TestPipelines;
}

my $output_dir = get_output_dir('breakdancer_analysis_pipeline');

ok my $pipeline = VRPipe::Pipeline->get(name => 'breakdancer_analysis'), 'able to get the breakdancer_analysis pipeline';
my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(breakdancer_bam2cfg breakdancer_sv_detection);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->get(name => 'my breakdancer_analysis pipeline setup',
		datasource => VRPipe::DataSource->get(type => 'delimited',
			method => 'all_columns',
			options => { delimiter => "\t" },
			source => file(qw(t data hs_chr20.bam.datasource))),
		output_root => $output_dir,
		pipeline => $pipeline,
		options => { 
		    'bam2cfg_options' => '-q 20 -c 3 -n 100000',
		    'breakdancer_max_options' => '-m 10000000 -q 25 -y 20',
			'whole_genome_mode' => 1,
			#'chrom_list' => '20',
		    cleanup => 0,
        });

my (@output_files,@final_files);
my $element_id=0;
foreach my $in ('hs_chr20.a', 'hs_chr20.c') {
    $element_id++;
    my @output_dirs = output_subdirs($element_id);
    push(@output_files, file(@output_dirs, '1_breakdancer_bam2cfg', "bd.cfg"));
    push(@final_files, file(@output_dirs, '2_breakdancer_sv_detection', "bd.max"));
}
ok handle_pipeline(@output_files,@final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
