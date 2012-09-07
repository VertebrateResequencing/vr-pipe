#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES GATK)]);
    use TestPipelines;
    use_ok('VRPipe::Steps::gatk_left_align_variants');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('gatk_left_align_variants', 'vcf_files');
is_deeply [$step->id, $step->description], [1, 'Runs GATK LeftAlignVariants to left-align indels in VCF files'], 'gatk_left_align_variants step created and has correct description';

my $original_ref_fa = VRPipe::File->create(path => file(qw(t data human_g1k_v37.chr11.chr20.fa.gz))->absolute);
my $ref_fa = VRPipe::File->create(path => file($output_dir, 'human_g1k_v37.chr11.chr20.fa')->absolute);
my $oh = $original_ref_fa->openr;
my $nh = $ref_fa->openw;
while (<$oh>) {
    print $nh $_;
}
close($oh);
close($nh);

my $fofn = VRPipe::File->create(path => file($output_dir, 'left_align.fofn'));
my $fh = $fofn->openw;
print $fh file(qw(t data dbsnp_132.b37.chr20_reduced.vcf.gz))->absolute->stringify . "\n";
$fofn->close;

my $setup = VRPipe::PipelineSetup->create(
    name       => 'gatk_left_align_variants',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn',
        method  => 'all',
        options => {},
        source  => $fofn->path
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta             => $ref_fa->path,
        left_align_variants_options => '-L 20:1-600000'
    }
);

my @output_subdirs = output_subdirs(1);
ok handle_pipeline(file(@output_subdirs, '1_gatk_left_align_variants', "dbsnp_132.b37.chr20_reduced.aln.0.vcf.gz")), 'gatk_left_align_variants pipeline ran ok, generating the expected file';
