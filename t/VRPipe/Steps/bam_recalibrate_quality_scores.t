#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES GATK)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::bam_recalibrate_quality_scores');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_recalibrate_quality_scores', 'bam_files');
VRPipe::StepAdaptor->create(pipeline => $pipeline, to_step => 1, adaptor_hash => { 'bai_files' => { data_element => 0 }, 'bam_files' => { data_element => 0 }, 'bam_recalibration_files' => { data_element => 0 } });
is_deeply [$step->id, $step->description], [1, 'Recalibrate quality scores of each mapped base using GATK'], 'bam_recalibrate_quality_scores step created and has correct description';

my $ref_fa = file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify;
my $bam    = VRPipe::File->create(path => file(qw(t data SRR790038.se.realign.bam))->absolute->stringify, type => 'bam');
my $bai    = VRPipe::File->create(path => file(qw(t data SRR790038.se.realign.bam.bai))->absolute->stringify, type => 'bin');
my $recal  = VRPipe::File->create(path => file(qw(t data SRR790038.se.realign.recal_data.csv))->absolute->stringify, type => 'txt');
$recal->add_metadata({ source_bam => file(qw(t data SRR790038.se.realign.bam))->absolute->stringify });

# test as part of a pipeline
my $setup = VRPipe::PipelineSetup->create(
    name        => 'recal_step_test',
    datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'group_all', source => file(qw(t data recal_with_few_reads.fofn))->absolute),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options => { reference_fasta => $ref_fa, bam_recalibration_options => '-l INFO --disable_bam_indexing' }
);

my @output_subdirs = output_subdirs(1);
ok handle_pipeline(file(@output_subdirs, '1_bam_recalibrate_quality_scores', qq[SRR790038.se.realign.recal.bam])), 'bam_recalibrate_quality_scores single-step pipeline ran ok and created expected output files';

finish;
