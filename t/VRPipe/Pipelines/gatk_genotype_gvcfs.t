#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES GATK GATK_JAVA)],
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('gatk_genotype_gvcfs');

# check pipeline has correct steps
ok my $gvcf_merge_pipeline = VRPipe::Pipeline->create(name => 'gatk_genotype_gvcfs'), 'able to get the gatk_genotype_gvcfs pipeline';
my @sb_names;
foreach my $stepmember ($gvcf_merge_pipeline->step_members) {
    push(@sb_names, $stepmember->step->name);
}

is_deeply \@sb_names, [qw(gatk_combine_gvcfs gatk_genotype_gvcfs)], 'the gatk_genotype_gvcfs pipeline has the correct steps';

my $merge_setup = VRPipe::PipelineSetup->create(
    name       => 'GATK merge gVCFs setup',
    pipeline   => $gvcf_merge_pipeline,
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'group_all',
        source => file(qw(t data gvcfs_to_merge.fofn)),
    ),
    options => {
        reference_fasta          => file(qw(t data GRCh38.chr1.fna))->absolute,
        maximum_gvcfs_to_combine => 3,
        gatk_path                => $ENV{GATK},
        java_exe                 => $ENV{GATK_JAVA},
    },
    output_root => $output_dir
);

my @final_files;
foreach my $element (@{ get_elements($merge_setup->datasource) }) {
    my @output_dirs = output_subdirs($element->id, $merge_setup->id);
    push(@final_files, file(@output_dirs, '2_gatk_genotype_gvcfs', 'gatk_mergeGvcf.vcf.gz'));
    push(@final_files, file(@output_dirs, '2_gatk_genotype_gvcfs', 'gatk_mergeGvcf.vcf.gz.tbi'));
}
ok handle_pipeline(@final_files), 'gatk_genotype_gvcfs pipeline ran ok';

finish;
exit;
