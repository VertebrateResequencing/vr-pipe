#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES GATK)]);
    use TestPipelines;
    use_ok('VRPipe::Steps::gatk_variant_annotator');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('gatk_variant_annotator', 'vcf_files');
is_deeply [$step->id, $step->description], [1, 'Run GATK VariantAnnotator to annotate multiple VCF files'], 'gatk_variant_annotator step created and has correct description';

my $original_ref_fa = VRPipe::File->create(path => file(qw(t data human_g1k_v37.chr11.chr20.fa.gz))->absolute);
my $ref_fa = VRPipe::File->create(path => file($output_dir, 'human_g1k_v37.chr11.chr20.fa')->absolute);
my $oh = $original_ref_fa->openr;
my $nh = $ref_fa->openw;
while (<$oh>) {
    print $nh $_;
}
close($oh);
close($nh);

my $fofn = VRPipe::File->create(path => file($output_dir, 'bams.list'));
my $fh = $fofn->openw;
print $fh file(qw(t data chrom20.ILLUMINA.bwa.CEU.low_coverage.bam))->absolute->stringify . "\n";
print $fh file(qw(t data chrom20.ILLUMINA.bwa.JPT.low_coverage.bam))->absolute->stringify . "\n";
print $fh file(qw(t data chrom20.SOLID.bfast.JPT.low_coverage.bam))->absolute->stringify . "\n";
close($fh);

my $setup = VRPipe::PipelineSetup->create(
    name       => 'gatk_variant_annotator',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn',
        method  => 'all',
        options => {},
        source  => file(qw(t data annotation.fofn))->absolute
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta           => $ref_fa->path,
        variant_annotator_options => '-L 20:1-75500 --group StandardAnnotation -I ' . $fofn->path
    }
);
my @output_subdirs = output_subdirs(1);
ok handle_pipeline(file(@output_subdirs, '1_gatk_variant_annotator', "annotation.anno.vcf.gz")), 'gatk_variant_annotator pipeline ran ok, generating the expected file';
