#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => 'VRPIPE_TEST_PIPELINES',
		    required_exe => 'glf');
    use TestPipelines;
}

my $checking_output_dir = get_output_dir('genotype_checking_wgs_pipeline');

ok my $gc_pipeline = VRPipe::Pipeline->get(name => 'genotype_checking_wgs'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($gc_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(snp_bin_hapmap_sites
                             mpileup_bcf_hapmap
							 glf_check_genotype
							 gtypex_genotype_analysis);

is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $snp_bin_source = file(qw(t data qc_chr20_snps.bin));
my $snp_dir = dir($checking_output_dir, 'snp');
$gc_pipeline->make_path($snp_dir);
my $snp_bin = file($snp_dir, 'qc_chr20_snps.bin')->stringify;
copy($snp_bin_source, $snp_bin);

my $ref_fa_source = file(qw(t data human_g1k_v37.chr20.fa));
my $ref_dir = dir($checking_output_dir, 'ref');
$gc_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'human_g1k_v37.chr20.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $mpileup_opts = '-ugDI -d 1000';

my $gc_pipelinesetup = VRPipe::PipelineSetup->get(name => 'genotype_checking',
                                                  datasource => VRPipe::DataSource->get(type => 'fofn',
                                                                                        method => 'all',
                                                                                        source => file(qw(t data hs_chr20.qc.bam.fofn))->absolute),
                                                  output_root => $checking_output_dir,
                                                  pipeline => $gc_pipeline,
                                                  options => {reference_fasta => $ref_fa,
                                                              snp_binary_file => $snp_bin,
                                                              samtools_mpileup_options => $mpileup_opts});

my (@output_files,@final_files);
my $element_id=0;
foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@output_files, file(@output_subdirs, '2_mpileup_bcf_hapmap', "${in}.bam.bcf"));
    push(@final_files, file(@output_subdirs, '3_glf_check_genotype', "${in}.bam.bcf.gtypex"));
    push(@final_files, file(@output_subdirs, '4_gtypex_genotype_analysis', "${in}.bam.gtype"));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;

finish;