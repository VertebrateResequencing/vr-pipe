#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(samtools bcftools htscmd)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('bam_htscmd_genotype_checking');

ok my $pipeline = VRPipe::Pipeline->create(name => 'bam_htscmd_genotype_checking'), 'able to get the bam_htscmd_genotype_checking pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(bam_index vcf_sites mpileup_vcf vcf_index htscmd_gtcheck htscmd_genotype_analysis);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $vcf_geno_source = file(qw(t data hs_chr20.genotypes.vcf.gz));
my $vcfg_dir = dir($output_dir, 'vcf_geno');
$pipeline->make_path($vcfg_dir);
my $vcf_geno = file($vcfg_dir, 'hs_chr20.genotypes.vcf.gz')->stringify;
copy($vcf_geno_source,       $vcf_geno);
copy("$vcf_geno_source.tbi", "$vcf_geno.tbi");

my $ref_fa_source = file(qw(t data human_g1k_v37.chr20.fa));
my $ref_dir = dir($output_dir, 'ref');
$pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'human_g1k_v37.chr20.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $ds = VRPipe::DataSource->create(
    type   => 'fofn',
    method => 'all',
    source => file(qw(t data hs_chr20.qc.bam.fofn))->absolute
);

VRPipe::PipelineSetup->create(
    name        => 'my bam_htscmd_genotype_checking pipeline setup',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta => $ref_fa,
        genotypes_vcf   => $vcf_geno,
        cleanup         => 0
    }
);

my (@output_files, @final_files);
my $element_id = 0;
foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
    # 3 of the input samples are actually NA20544 while 1 is NA20588; for
    # testing purposes we set them all to NA20544 so that we get a failure
    # (the SM tag in the bam headers is wrong?)
    VRPipe::File->create(path => file('t', 'data', $in . '.bam')->absolute)->add_metadata({ sample => 'NA20544' });
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@output_files, file(@output_subdirs, '3_mpileup_vcf',    "mpileup.vcf.gz"));
    push(@final_files,  file(@output_subdirs, '5_htscmd_gtcheck', "mpileup.vcf.gz.gtypex"));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

my @found_gtype_analysis;
foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
    my $meta = VRPipe::File->create(path => file('t', 'data', $in . '.bam')->absolute)->metadata;
    my ($s, $e, $f, $r, $c) = $meta->{gtype_analysis} =~ /status=(\S+) expected=(\S+) found=(\S+) ratio=(\S+) concordance=(\S+)/;
    my $r_ok = $r >= 1  ? 1 : 0;
    my $c_ok = $c > 0.3 ? 1 : 0;
    push(@found_gtype_analysis, [$s, $e, $f, $r_ok, $c_ok]);
}
my @expected_gtype_analysis = (
    [qw(wrong NA20544 NA20588 1 0)],
    [qw(confirmed NA20544 NA20544 1 1)],
    [qw(confirmed NA20544 NA20544 1 1)],
    [qw(confirmed NA20544 NA20544 1 1)],
);

is_deeply \@found_gtype_analysis, \@expected_gtype_analysis, 'pipeline generated the expected genotype analysis metadata';

# also test the multiple_samples_per_individual mode
$output_dir = get_output_dir('bam_htscmd_genotype_checking_multiple_samples_per_individual');
VRPipe::PipelineSetup->create(
    name        => 'my bam_htscmd_genotype_checking multiple_samples_per_individual pipeline setup',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta                 => $ref_fa,
        genotypes_vcf                   => $vcf_geno,
        min_concordance                 => 0.4,
        multiple_samples_per_individual => 1,
        cleanup                         => 0
    }
);

(@output_files, @final_files) = ();
$element_id = 0;
foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 2);
    push(@output_files, file(@output_subdirs, '3_mpileup_vcf',    "mpileup.vcf.gz"));
    push(@final_files,  file(@output_subdirs, '5_htscmd_gtcheck', "mpileup.vcf.gz.gtypex"));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran in multiple_samples_per_individual mode and created all expected output files';

@found_gtype_analysis = ();
foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
    my $meta = VRPipe::File->create(path => file('t', 'data', $in . '.bam')->absolute)->metadata;
    my ($s, $e, $f, $c, $r) = $meta->{gtype_analysis} =~ /status=(\S+) expected=(\S+) found=(\S+) concordance=(\S+) ratio=(\S+)/;
    my $c_ok = $c > 0.3 ? 1 : 0;
    my $r_ok = $r >= 1  ? 1 : 0;
    push(@found_gtype_analysis, [$s, $e, $f, $c_ok, $r_ok]);

}
@expected_gtype_analysis = (
    [qw(wrong NA20544 NA20588 0 1)],
    [qw(confirmed NA20544 NA20544 1 1)],
    [qw(confirmed NA20544 NA20544 1 1)],
    [qw(confirmed NA20544 NA20544 1 1)],
);

is_deeply \@found_gtype_analysis, \@expected_gtype_analysis, 'pipeline generated the expected genotype analysis metadata';

done_testing;
exit;
