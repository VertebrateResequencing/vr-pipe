#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => 'verifyBamID'
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('verify_bamid_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'verify_bamid'), 'able to get the verify_bamid pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(verify_bamid);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $fofn = file(qw(t data bams.fofn))->absolute->stringify;
my $vcf  = file(qw(t data verify_ref.vcf.gz))->absolute->stringify;

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my verify_bamid pipeline setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => $fofn
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        verify_bamid_opts => "--vcf $vcf --ignoreRG",
        cleanup           => 0
    }
);

my (@output_files);
my $element_id = 0;
foreach my $sample (qw(NA19334 NA19381 NA20281)) {
    $element_id++;
    my @output_dirs = output_subdirs($element_id);
    push(@output_files, file(@output_dirs, '1_verify_bamid', "$sample.log"));
    push(@output_files, file(@output_dirs, '1_verify_bamid', "$sample.selfSM"));
    push(@output_files, file(@output_dirs, '1_verify_bamid', "$sample.depthSM"));
}

ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
