#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use_ok('VRPipe::Persistent::Schema');
    use TestPipelines;
}

my $output_dir = get_output_dir('snp_calling_gatk_vcf_pipeline');

ok my $pipeline = VRPipe::Pipeline->get(name => 'snp_calling_gatk_vcf'), 'able to get the snp_calling_gatk_vcf pipeline';

my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(chunk_genomic_region gatk_genotype vcf_concat);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->get(name => 'my snp_calling_gatk_vcf pipeline setup',
		datasource => VRPipe::DataSource->get(type => 'delimited',
			method => 'all_columns',
			options => { delimiter => "\t" },
			source => file(qw(t data hs_chr20.bam.fofn))),
		output_root => $output_dir,
		pipeline => $pipeline,
		options => { cleanup => 0, 
			#chunking_regions_file => file(qw(t data human_g1k_v37.fasta.fai))->absolute->stringify,
			chunking_regions_file => file(qw(t data hs_chunking_regions.list))->absolute->stringify,
			dbsnp_ref => file(qw(t data dbsnp_130_b37.chr20.reduced.rod))->absolute->stringify,
			reference_fasta => file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify,
			interval_list => file(qw(t data hs_chr20.invervals.list))->absolute->stringify,
			genotyper_opts => '-stand_call_conf 50.0 -stand_emit_conf 10.0',
			#genotyper_opts => '-stand_call_conf 50.0 -stand_emit_conf 10.0 -L 20:1-70000', 
		}
	);

my (@output_files,@final_files);
my @files = ('hs_chr20.a.bam','hs_chr20.c.bam');
my $element_id = 0;
foreach (@files) {
  $element_id++;
  my $file = 'merged.vcf.gz';
  push(@output_files, file($output_dir, output_subdirs($element_id), '3_vcf_concat', $file));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
