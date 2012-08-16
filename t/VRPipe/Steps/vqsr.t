#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES GATK)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::gatk_variant_recalibration');
    use_ok('VRPipe::Steps::gatk_variant_recalibration_for_snps');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('gatk_variant_recalibration', 'vcf_files');
is_deeply [$step->id, $step->description], [1, 'Recalibrates variant calls using Variant Quality Score Recalibration (VQSR)'], 'gatk_variant_recalibration step created and has correct description';

my $recal_opts = "--percentBadVariants 0.07 --maxGaussians 1 --minNumBadVariants 10";
$recal_opts .= " -resource:dbsnp,known=true,training=true,truth=false,prior=8.0 " . file(qw(t data dbsnp_132.b37.chr20_reduced.vcf.gz))->absolute->stringify;
$recal_opts .= " -resource:omni,known=false,training=true,truth=false,prior=12.0 " . file(qw(t data 1000G_omni2.5.b37.sites.chr20_reduced.vcf.gz))->absolute->stringify;
$recal_opts .= " -resource:hapmap,known=false,training=true,truth=true,prior=15.0 " . file(qw(t data hapmap_3.3.b37.sites.chr20_reduced.vcf.gz))->absolute->stringify;
$recal_opts .= ' -an QD -an HaplotypeScore';

VRPipe::PipelineSetup->create(
    name        => 'vqsr step test',
    datasource  => VRPipe::DataSource->create(type => 'fofn_with_metadata', method => 'group_all', source => file(qw(t data vqsr_datasource.fofn))->absolute),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta               => file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify,
        variant_recalibration_options => $recal_opts
    }
);

my @output_files;
my @output_subdirs = output_subdirs(1, 1);
foreach my $file (qw(recal recal.tranches recal.tranches.pdf recal.r recal.r.pdf)) {
    push(@output_files, file(@output_subdirs, '1_gatk_variant_recalibration', 'BOTH.' . $file));
}

ok handle_pipeline(@output_files), 'single step gatk_variant_recalibration pipeline ran and created all expected output files';

($output_dir, $pipeline, $step) = create_single_step_pipeline('gatk_variant_recalibration_for_snps', 'vcf_files');
is_deeply [$step->id, $step->description], [2, 'Recalibrates SNP calls using Variant Quality Score Recalibration (VQSR)'], 'gatk_variant_recalibration_for_snps step created and has correct description';

VRPipe::PipelineSetup->create(
    name        => 'vqsr step test',
    datasource  => VRPipe::DataSource->create(type => 'fofn_with_metadata', method => 'group_all', source => file(qw(t data vqsr_datasource.fofn))->absolute),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta           => file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify,
        snp_recalibration_options => $recal_opts
    }
);
@output_files = ();
@output_subdirs = output_subdirs(1, 2);
foreach my $file (qw(recal recal.tranches recal.tranches.pdf recal.r recal.r.pdf)) {
    push(@output_files, file(@output_subdirs, '1_gatk_variant_recalibration_for_snps', 'SNP.' . $file));
}
ok handle_pipeline(@output_files), 'single step gatk_variant_recalibration_for_snps pipeline ran and created all expected output files';

finish;
