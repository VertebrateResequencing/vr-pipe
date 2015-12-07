#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES GATK GATK2 JAVA7)],
        required_exe => [qw(java samtools bammarkduplicates)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('bam_improvement_gatk_v1_merge_and_mark_duplicates');
ok my $pipeline = VRPipe::Pipeline->create(name => 'bam_improvement_gatk_v1_merge_and_mark_duplicates'), 'able to get the bam_improvement_gatk_v1_merge_and_mark_duplicates pipeline';

my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(samtools_add_readgroup gatk_target_interval_creator bam_realignment_around_known_indels bam_count_covariates bam_recalibrate_quality_scores samtools_merge_and_bam_mark_duplicates);
is_deeply \@s_names, \@expected_step_names, 'the pipelines have the correct steps';

my $res_dir = dir($output_dir, 'resources');
$pipeline->make_path($res_dir);

my $known_indels_source = file(qw(t data known_indels_pombe.vcf.gz));
my $known_indels = file($res_dir, 'known_indels.vcf.gz')->stringify;
copy($known_indels_source,          $known_indels);
copy($known_indels_source . '.tbi', $known_indels . '.tbi');

my $known_sites_source = file(qw(t data known_sites_pombe.vcf.gz));
my $known_sites = file($res_dir, 'known_sites.vcf.gz')->stringify;
copy($known_sites_source,          $known_sites);
copy($known_sites_source . '.tbi', $known_sites . '.tbi');

my $fofn_file = VRPipe::File->create(path => file(qw(t data pombe_bams.fofn))->absolute);

VRPipe::PipelineSetup->create(
    name       => 'bam improvement and merge',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => $fofn_file->path->stringify,
        options => { metadata_keys => 'library' }
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        java_exe                      => 'java',
        samtools_exe                  => 'samtools',
        gatk_path                     => $ENV{GATK},
        reference_fasta               => file(qw(t data pombe_ref.fa))->absolute->stringify,
        known_indels_for_realignment  => "-known $known_indels",
        known_sites_for_recalibration => "-knownSites $known_sites",
        gatk_count_covariates_options => '-l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate',
        bam_realignment_options       => '-LOD 0.4 -model KNOWNS_ONLY -compress 0',
        bam_recalibration_options     => '-l INFO',
        cleanup                       => 0,
    }
);

ok handle_pipeline(), 'bam_improvement_gatk_v1_merge_and_mark_duplicates pipeline ran okay';

## test the same pipeline based on GATK v2.x
$output_dir = get_output_dir('bam_improvement_gatk_v2_merge_and_mark_duplicates');
ok $pipeline = VRPipe::Pipeline->create(name => 'bam_improvement_gatk_v2_merge_and_mark_duplicates'), 'able to get the bam_improvement_gatk_v2_merge_and_mark_duplicates pipeline';

@s_names = ();
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
@expected_step_names = qw(samtools_add_readgroup gatk_target_interval_creator bam_realignment_around_known_indels gatk_base_recalibrator gatk_print_reads_with_bqsr samtools_merge_and_bam_mark_duplicates);
is_deeply \@s_names, \@expected_step_names, 'the pipelines have the correct steps';

VRPipe::PipelineSetup->create(
    name       => 'bam improvement and merge GATK v2',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => $fofn_file->path->stringify,
        options => { metadata_keys => 'library' }
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        java_exe                     => $ENV{JAVA7},
        samtools_exe                 => 'samtools',
        gatk_path                    => $ENV{GATK2},
        reference_fasta              => file(qw(t data pombe_ref.fa))->absolute->stringify,
        known_indels_for_realignment => "-known $known_indels",
        bam_realignment_options      => '-LOD 0.4 -model KNOWNS_ONLY -compress 0',
        base_recalibrator_options    => "-l INFO -knownSites $known_sites -cov CycleCovariate -cov ContextCovariate",
        print_reads_options          => '-l INFO',
        cleanup                      => 0,
    }
);

ok handle_pipeline(), 'bam_improvement_gatk_v2_merge_and_mark_duplicates pipeline ran okay';

finish;
