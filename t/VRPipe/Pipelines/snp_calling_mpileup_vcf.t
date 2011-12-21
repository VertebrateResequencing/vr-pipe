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
		datasource => VRPipe::DataSource->get(type => 'delimited',
			method => 'all_columns',
			options => { delimiter => "\t" },
			source => file(qw(t data hs_chr20.bam.fofn))),
		output_root => $output_dir,
		pipeline => $pipeline,
		options => { cleanup => 0,
			ploidy_definition => "{default=>2,X=>[{region=>'1-60000',M=>1},{region=>'2699521-154931043',M=>1},],Y=>[{region=>'1-59373566',M=>1,F=>0},],}",
			#interval_list => file(qw(t data hs_chr20.invervals.bed))->absolute->stringify,
			#samtools_mpileup_options => '-C50 -aug -r 20:1-70000',
			samtools_mpileup_options => '-C50 -aug',
			reference_fasta => file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify,
			#reference_fasta => "/lustre/scratch105/projects/g1k/ref/main_project/human_g1k_v37.fasta",
		}
);


my (@output_files,@final_files);
my @files = ('hs_chr20.a.bam','hs_chr20.c.bam');
my $element_id = 0;
foreach (@files) {
  $element_id++;
  my $file = 'mpileup.vcf.gz';
  push(@output_files, file($output_dir, output_subdirs($element_id), '1_mpileup_vcf', $file));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
