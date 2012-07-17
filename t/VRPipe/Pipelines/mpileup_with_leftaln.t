#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],
		    required_exe => [qw(samtools vcf-concat bcftools)]);
    use TestPipelines;
}

my $output_dir = get_output_dir('mpileup_with_leftaln_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'mpileup_with_leftaln'), 'able to get the mpileup_with_leftaln pipeline';
my @s_names;

foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(mpileup_vcf vcf_index gatk_vcf_leftalign);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->create(name => 'my mpileup_with_leftaln pipeline setup',
		datasource => VRPipe::DataSource->create(type => 'fofn',
			method => 'all',
			source => file(qw(t data hs_chr20.bam.fofn))),
		output_root => $output_dir,
		pipeline => $pipeline,
		options => { cleanup => 0,
			samtools_mpileup_options => '-C50 -aug',
			reference_fasta => file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify,
		}
);

my (@output_files,@final_files);
my $element_id=0;
foreach my $in ('hs_chr20.a', 'hs_chr20.b') {
    $element_id++;
    my @output_dirs = output_subdirs($element_id);
    push(@final_files, file(@output_dirs, '3_gatk_vcf_leftalign', "${in}.mpileup.aln.vcf.gz"));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
