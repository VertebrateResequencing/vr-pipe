#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Data::Dumper;

BEGIN {
    use Test::Most tests => 4;
    use_ok('VRPipe::Persistent::Schema');
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
		}
);


my (@output_files,@final_files);
my @files = ('hs_chr20.a.bam','hs_chr20.c.bam');
my $element_id = 0;
foreach (@files) {
  $element_id++;
  my $file = 'mpileup.vcf';
  push(@output_files, file($output_dir, output_subdirs($element_id), 'mpileup_vcf', $file));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
