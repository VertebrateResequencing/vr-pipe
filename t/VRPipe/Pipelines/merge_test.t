#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES PICARD)],
        required_exe => [qw(samtools bamcheck)]
    );
    use TestPipelines;
}

my $merge_output_dir = get_output_dir('merge_test');

ok my $merge_pipeline = VRPipe::Pipeline->create(name => 'merge_test_pipeline'), 'able to get the merge_test_pipeline pipeline';

my @s_names;
foreach my $stepmember ($merge_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(test_import_bams
  bam_metadata
  bam_strip_tags
  bam_merge
  bam_mark_duplicates);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

ok my $ds = VRPipe::DataSource->create(
    type    => 'delimited',
    method  => 'grouped_single_column',
    source  => file(qw(t data datasource.bams))->absolute->stringify,
    options => {
        delimiter => "\t",
        group_by  => 1,
        column    => 2
    }
  ),
  'could create a delimited datasource';

my @results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply \@results, [{ paths => [file('t', 'data', '2822_6.pe.bam')->absolute, file('t', 'data', '2822_6.se.bam')->absolute, file('t', 'data', '2822_7.pe.bam')->absolute], group => 'LIB01' }, { paths => [file('t', 'data', '2823_4.pe.bam')->absolute], group => 'LIB02' }, { paths => [file('t', 'data', '8324_8.pe.bam')->absolute], group => 'LIB03' }], 'got correct results for delimited datasource';

my $merge_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 's_suis merge',
    datasource => VRPipe::DataSource->create(
        type    => 'delimited',
        method  => 'grouped_single_column',
        source  => file(qw(t data datasource.bams))->absolute->stringify,
        options => {
            delimiter => "\t",
            group_by  => 1,
            column    => 2
        }
    ),
    output_root => $merge_output_dir,
    pipeline    => $merge_pipeline,
    options     => {
        reference_fasta                       => file(qw(t data S_suis_P17.fa))->absolute->stringify,
        bam_tags_to_strip                     => 'OQ XM XG XO',
        bam_merge_keep_single_paired_separate => 1,
        cleanup                               => 0
    }
);

ok handle_pipeline(), 'pipeline ran ok';

is_deeply [VRPipe::StepState->create(pipelinesetup => 1, stepmember => 4, dataelement => 1)->cmd_summary->summary, VRPipe::StepState->create(pipelinesetup => 1, stepmember => 5, dataelement => 1)->cmd_summary->summary], ['java $jvm_args -jar MergeSamFiles.jar INPUT=$bam_file(s) OUTPUT=$merged_bam VALIDATION_STRINGENCY=SILENT', 'java $jvm_args -jar MarkDuplicates.jar INPUT=$bam_file OUTPUT=$markdup_bam_file ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT'], 'cmd summaries for the major steps were as expected';

finish;
