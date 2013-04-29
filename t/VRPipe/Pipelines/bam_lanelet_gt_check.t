#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('bam_lanelet_gt_check');

ok my $pipeline = VRPipe::Pipeline->create(name => 'bam_lanelet_gt_check'), 'able to get the bam_lanelet_gt_check pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(lanelet_gt_bam_select bam_index vcf_sites mpileup_vcf vcf_index htscmd_gtcheck lanelet_gt_bam_update);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $vcf_geno_source = file(qw(t data hs_chr20.genotypes.vcf.gz));
my $vcfg_dir = dir($output_dir, 'vcf_geno');
$pipeline->make_path($vcfg_dir);
my $vcf_geno = file($vcfg_dir, 'hs_chr20.genotypes.vcf.gz')->stringify;
copy($vcf_geno_source, $vcf_geno);
copy("$vcf_geno_source.tbi", "$vcf_geno.tbi");

my $ref_fa_source = file(qw(t data human_g1k_v37.chr20.fa));
my $ref_dir = dir($output_dir, 'ref');
$pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'human_g1k_v37.chr20.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $ds = VRPipe::DataSource->create(
    type   => 'fofn_with_metadata',
    method => 'grouped_by_metadata',
    options => { metadata_keys => 'library'},
    source => file(qw(t data hs_chr20.gtcheck.bam.fofnwmd))
);

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name        => 'my bam_lanelet_gt_check pipeline setup',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        minimum_lib_gbm => 5,
        minimum_gbm_confirmed => 5,
        reference_fasta => $ref_fa,
        genotypes_vcf => $vcf_geno,
        cleanup => 0
    }
);

my (@output_files);
for (my $element_id=1; $element_id<3; $element_id++) {
    my @output_subdirs = output_subdirs($element_id);
    push(@output_files, file(@output_subdirs, '4_mpileup_vcf', "mpileup.vcf.gz"));
    push(@output_files, file(@output_subdirs, '6_htscmd_gtcheck', "mpileup.vcf.gz.gtypex"));
}

ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

my @bam_meta;
foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
    my $bam = VRPipe::File->get(path => file('t', 'data', "$in.bam")->absolute);
    push (@bam_meta,$bam->metadata->{library_gt});

}
my @expected_meta=(undef,undef,"1",undef);
is_deeply \@expected_meta, \@bam_meta, 'the bam metadata library_gt was updated correctly';

done_testing;
exit;
