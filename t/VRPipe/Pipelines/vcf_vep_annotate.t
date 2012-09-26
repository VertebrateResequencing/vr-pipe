#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(variant_effect_predictor.pl vcf2consequences_vep tabix)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('vcf_vep_annotate_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'vcf_vep_annotate'), 'able to get the vcf_vep_annotate pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(vep_analysis vcf_vep_consequences vcf_index);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $vep_cache  = file(qw(t data vep_cache))->absolute->stringify;
my $gerp_cache = file(qw(t data gerp_cache))->absolute->stringify;

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my vcf_vep_annotate pipeline setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data datasource.vcf_fofn))
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        'vep_options'              => "--sift b --polyphen b --condel b --gene --hgnc --format vcf --force_overwrite --offline --cache --dir $vep_cache",
        'vcf2consequences_options' => "--gerp $gerp_cache",
        cleanup                    => 1
    }
);

my (@output_files, @final_files);
my $element_id = 0;
foreach my $in ('test1', 'test2') {
    $element_id++;
    my @output_dirs = output_subdirs($element_id);
    push(@final_files, file(@output_dirs, '2_vcf_vep_consequences', "${in}.conseq.vcf.gz"));
    push(@final_files, file(@output_dirs, '2_vcf_vep_consequences', "${in}.conseq.vcf.gz.tbi"));
}
ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
