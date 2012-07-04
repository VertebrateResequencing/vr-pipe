#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)], );
    use TestPipelines;
}

my $output_dir = get_output_dir('vcf_chunked_vep_annotate_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'vcf_chunked_vep_annotate'), 'able to get the vcf_chunked_vep_annotate pipeline';
my @s_names;

foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(chunk_genomic_region vcf_split vep_analysis vcf_vep_consequences vcf_concat vcf_index);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $vep_cache = file(qw(t data vep_cache))->absolute->stringify;
my $gerp_cache = file(qw(t data gerp_cache))->absolute->stringify;

my $test_pipelinesetup = VRPipe::PipelineSetup->create(name => 'my vcf_chunked_vep_annotate pipeline setup',
		datasource => VRPipe::DataSource->create(type => 'fofn',
			method => 'all',
			source => file(qw(t data datasource.vcf_fofn))),
		output_root => $output_dir,
		pipeline => $pipeline,
		options => { cleanup => 0,
			genomic_region_file => file(qw(t data human_g1k_v37.fasta.fai))->absolute->stringify,
			chunk_size => 20000000,
			'vep_options' => "--sift b --polyphen b --condel b --gene --hgnc --format vcf --force_overwrite --cache --dir $vep_cache",
			'vcf2consequences_options' => "--grantham --gerp $gerp_cache",
		}
);

my (@output_files,@final_files);
my @files = ('test1.vcf.gz','test2.vcf.gz ');
my $element_id = 0;
foreach (@files) {
  $element_id++;
  my @output_subdirs = output_subdirs($element_id);
  my $file = 'merged.vcf.gz';
  push(@output_files, file(@output_subdirs, '5_vcf_concat', $file));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
