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

my $output_dir = get_output_dir('snp_calling_mpileup_vcf_pipeline');

ok my $pipeline = VRPipe::Pipeline->get(name => 'snp_calling_mpileup_vcf'), 'able to get the snp_calling_mpileup_vcf pipeline';
my @s_names;

foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(mpileup_vcf);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->get(name => 'my snp_calling_mpileup_vcf pipeline setup',
		datasource => VRPipe::DataSource->get(type => 'fofn',
			method => 'all',
			source => file(qw(t data hs_chr20.bam.fofn))),
		output_root => $output_dir,
		pipeline => $pipeline,
		options => { cleanup => 0,
			#interval_list => file(qw(t data hs_chr20.invervals.bed))->absolute->stringify,
			#samtools_mpileup_options => '-C50 -aug -r 20:1-70000',
			samtools_mpileup_options => '-C50 -aug',
			reference_fasta => file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify,
			mimimum_calls => 0,
		}
);


my (@output_files,@final_files);
my @files = ('hs_chr20.a.bam','hs_chr20.b.bam');	# output file name and metadat based on last column of delimited data source
my $element_id = 0;
foreach my $f (@files) {
  $f =~ s/\.bam/.mpileup.vcf.gz/;
  $element_id++;
  push(@output_files, file(output_subdirs($element_id), '1_mpileup_vcf', "$f"));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
