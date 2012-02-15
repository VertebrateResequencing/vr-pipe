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

ok my $pipeline = VRPipe::Pipeline->get(name => 'vcf_chunked_vep_annotate'), 'able to get the vcf_chunked_vep_annotate pipeline';
my @s_names;

foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(chunk_genomic_region vcf_split vcf_annotate vep_analysis vcf_vep_consequences vcf_concat vcf_index);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $annot_file = file(qw(t data g1k_dbsnp132_annot.tab.gz))->absolute->stringify;
my $annot_desc_file = file(qw(t data g1k_dbsnp132_annot_desc.txt))->absolute->stringify;
my $annot_2_file = file(qw(t data annots-rsIDs-AFs.2011-10-05.tab.gz))->absolute->stringify;
my $annot_2_desc_file = file(qw(t data annots-rsIDs-AFs.2011-10-05.tab.gz.desc))->absolute->stringify;
my $vep_cache = file(qw(t data vep_cache))->absolute->stringify;
my $gerp_cache = file(qw(t data gerp_cache))->absolute->stringify;

my $test_pipelinesetup = VRPipe::PipelineSetup->get(name => 'my vcf_chunked_vep_annotate pipeline setup',
		datasource => VRPipe::DataSource->get(type => 'fofn',
			method => 'all',
			source => file(qw(t data datasource.vcf_fofn))),
		output_root => $output_dir,
		pipeline => $pipeline,
		options => { cleanup => 0,
			chunking_regions_file => file(qw(t data human_g1k_v37.fasta.fai))->absolute->stringify,
			chunk_size => 20000000,
			'vcf-annotate_options' => "-a $annot_file -d $annot_desc_file -c CHROM,FROM,REF,ALT,-,-,INFO/KGPilot123,INFO/dbSNP ",
			'vcf-annotate_2_options' => "-a $annot_2_file -d $annot_2_desc_file -c CHROM,POS,ID,REF,ALT,INFO/AF_AFR,INFO/AF_AMR,INFO/AF_ASN,INFO/AF_EUR,INFO/AF_MAX ",
			'vep_options' => "--sift b --polyphen b --condel b --gene --hgnc --format vcf --force_overwrite --cache --dir $vep_cache",
			'tmp_exe' => "/lustre/scratch106/user/cj5/vep2.2/vep_resort.sh",	# If bams split by chr, no re-sort required so this can just be 'cat'
			'vcf2consequences_options' => "-grantham --gerp $gerp_cache",
		}
);


my (@output_files,@final_files);
my @files = ('hs_chr20.a.bam','hs_chr20.c.bam');
my $element_id = 0;
foreach (@files) {
  $element_id++;
  my @output_subdirs = output_subdirs($element_id);
  my $file = 'merged.vcf.gz';
  push(@output_files, file(@output_subdirs, '6_vcf_concat', $file));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
