#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => 'VRPIPE_TEST_PIPELINES',
        required_exe => [qw(bwa samtools)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('bam_spatial_filter');

ok my $pipeline = VRPipe::Pipeline->create(name => 'bam_spatial_filter'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(calculate_bam_spatial_filter apply_bam_spatial_filter)], 'the pipeline has the correct steps';

# ref genome must be an absolute path
my $original_ref_fa = VRPipe::File->create(path => file(qw(t data S_suis_P17.fa))->absolute);
my $ref_fa = VRPipe::File->create(path => file($output_dir, 'S_suis_P17.fa')->absolute);
$original_ref_fa->copy($ref_fa);

my $pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'bam_spatial_filter',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn_with_metadata',
        method => 'grouped_by_metadata',
        options => {metadata_keys => 'lane'},
        source => file(qw(t data hs_bam.fofnwm))
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        spatial_filter_exe => '/software/solexa/bin/pb_calibration/v9.0/spatial_filter',
        reference_fasta  => $ref_fa->path,
        tag_number => 0,
        mark_qcfail => 0,
        cleanup => 0,
    }
);

my (@output_files);
my $element_id = 0;
foreach my $in ('2822_6.pe', '2822_7.pe') {
    $element_id++;
    my @output_dirs = output_subdirs($element_id);
    push(@output_files, file(@output_dirs, '1_calculate_bam_spatial_filter', "${in}.bam.filt"));
    push(@output_files, file(@output_dirs, '2_apply_bam_spatial_filter', "${in}.spfilt.bam"));
}
ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

finish;
