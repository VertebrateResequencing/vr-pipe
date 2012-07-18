#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],);
    use TestPipelines;
}

my $output_dir = get_output_dir('snp_calling_gatk_vcf_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'snp_calling_gatk_vcf'), 'able to get the snp_calling_gatk_vcf pipeline';

my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(gatk_genotype vcf_index gatk_variant_filter vcf_index gatk_recalibrate_variants gatk_apply_recalibration vcf_index);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

#my $recal_opts = "--mode 'BOTH' --ignore_filter 'HARD_TO_VALIDATE' --percentBadVariants 0.07 --maxGaussians 4 --target_titv '2.08'"; # production settings
my $recal_opts = "--mode 'BOTH' --ignore_filter 'HARD_TO_VALIDATE' --percentBadVariants 0.07 --maxGaussians 1 --minNumBadVariants 10"; # test with small no of variants
$recal_opts .= " -an QD -an HaplotypeScore";                                                                                           # in SNPS.pm $gatk->set_annotations('HaplotypeScore', 'MQRankSum', 'ReadPosRankSum', 'HRun');

$recal_opts .= " -resource:dbsnp,known=true,training=true,truth=false,prior=8.0 " . file(qw(t data dbsnp_132.b37.chr20_reduced.vcf.gz))->absolute->stringify;
$recal_opts .= " -resource:omni,known=false,training=true,truth=false,prior=12.0 " . file(qw(t data 1000G_omni2.5.b37.sites.chr20_reduced.vcf.gz))->absolute->stringify;
$recal_opts .= " -resource:hapmap,known=false,training=true,truth=true,prior=15.0 " . file(qw(t data hapmap_3.3.b37.sites.chr20_reduced.vcf.gz))->absolute->stringify;
$recal_opts .= " --phone_home NO_ET";

my $var_filter_opts = "--filterExpression 'MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1' --filterName HARD_TO_VALIDATE --mask:NAME,BED " . file(qw(t data chr20_trunc.pilot.indels.bed.mask))->absolute->stringify . " --maskName 'InDel' --clusterWindowSize 11 --phone_home NO_ET";

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my snp_calling_gatk_vcf pipeline setup',
    datasource => VRPipe::DataSource->create(type    => 'delimited',
                                             method  => 'all_columns',
                                             options => { delimiter => "\t" },
                                             source  => file(qw(t data hs_chr20.bam.datasource))),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        cleanup         => 0,
        dbsnp_ref       => file(qw(t data dbsnp_132.b37.chr20_reduced.vcf.gz))->absolute->stringify,
        reference_fasta => file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify,            # contig size tailored for test bams
        #interval_list => file(qw(t data hs_chr20.invervals.list))->absolute->stringify,
        #interval_list => file(qw(t data exome_B_GRCh37.intervals.window100_chr20.intervals))->absolute->stringify,
        genotyper_opts => "--genotype_likelihoods_model BOTH -stand_call_conf 10.0 --phone_home NO_ET", # -stand_call_conf 10.0 hardcoded into "GATK.pm"
        #genotyper_opts => '-stand_call_conf 50.0 -stand_emit_conf 10.0 -L 20:1-70000',
        var_filter_opts  => "$var_filter_opts",
        var_recal_opts   => "$recal_opts",
        apply_recal_opts => "--mode BOTH --ts_filter_level 100 --phone_home NO_ET", });

my (@output_files, @final_files);
my @files = ('hs_chr20.a.bam', 'hs_chr20.c.bam');
my $element_id = 0;
foreach (@files) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    my $file           = 'gatk_var.vcf.gz';
    push(@output_files, file(@output_subdirs, '1_gatk_genotype', $file));
    push(@output_files, file(@output_subdirs, '1_gatk_genotype', "$file.tbi"));
    $file = 'gatk_var.filt.vcf.gz';
    push(@output_files, file(@output_subdirs, '3_gatk_variant_filter', $file));
    push(@output_files, file(@output_subdirs, '3_gatk_variant_filter', "$file.tbi"));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
