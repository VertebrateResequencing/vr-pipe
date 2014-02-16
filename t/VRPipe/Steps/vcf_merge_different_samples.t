#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(bcftools)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::vcf_merge_different_samples');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('vcf_merge_different_samples', 'vcf_files');
is_deeply [$step->id, $step->description], [1, 'Merges compressed VCFs using bcftools merge which contain different samples to produce a single VCF containing all the input samples'], 'vcf_merge_different_samples step created and has correct description';

VRPipe::PipelineSetup->create(
    name        => 'vcf_merge_different_samples step test',
    datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'group_all', source => file(qw(t data vcf_merge_different_samples_datasource.fofn))->absolute),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);

my @output_files;
my @output_subdirs = output_subdirs(1, 1);
push(@output_files, file(@output_subdirs, '1_vcf_merge_different_samples', 'merged.vcf.gz'));

ok handle_pipeline(@output_files), 'single step vcf_merge_different_samples pipeline ran and created all expected output files';

my $vcf_file = shift @output_files;
open(my $fh, "zcat $vcf_file |");
my $samples;
while (<$fh>) {
    next unless /^#CHROM.+FORMAT\s+(.+)$/;
    $samples = $1;
    last;
}
close($fh);

is $samples, "SAMPLE01\tSAMPLE02", "merged VCF contains both input samples";

finish;
