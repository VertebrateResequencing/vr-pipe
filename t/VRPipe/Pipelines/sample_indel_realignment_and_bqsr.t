#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES GATK2 JAVA7)],
        required_exe => [qw(samtools)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('sample_indel_realignment_and_bqsr');
ok my $pipeline = VRPipe::Pipeline->create(name => 'sample_indel_realignment_and_bqsr'), 'able to get the sample_indel_realignment_and_bqsr pipeline';

my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(
  gatk_target_interval_creator_discovery
  bam_realignment_around_discovered_indels
  gatk_base_recalibrator
  gatk_print_reads_with_bqsr);
is_deeply \@s_names, \@expected_step_names, 'the pipelines have the correct steps';

my $res_dir = dir($output_dir, 'resources');
$pipeline->make_path($res_dir);

my $known_sites = file(qw(t data known_sites.vcf.gz))->absolute->stringify;

my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data improvement_datasource.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($output_dir, 'improvement_datasource.fofn'));
my $ifh = $orig_fofn_file->openr;
my $ofh = $fofn_file->openw;
my $idx = 0;
print $ofh "path\tindex\n";
while (<$ifh>) {
    chomp;
    print $ofh "$_\t$idx\n";
    print $ofh "$_.bai\t$idx\n";
    $idx++;
}
$orig_fofn_file->close;
$fofn_file->close;

my $bqsr_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 's_suis sample indel realignment and bqsr',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => $fofn_file->path->stringify,
        options => { metadata_keys => 'index' }
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        java_exe                    => $ENV{JAVA7},
        gatk_path                   => $ENV{GATK2},
        reference_fasta             => file(qw(t data S_suis_P17.fa))->absolute->stringify,
        gatk_indelrealigner_options => '-l INFO',
        base_recalibrator_options   => "-l INFO -knownSites $known_sites -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate",
        print_reads_options         => '-l INFO',
        cleanup                     => 0,
    }
);

ok handle_pipeline(), 'sample_indel_realignment_and_bqsr pipeline ran okay';

finish;
