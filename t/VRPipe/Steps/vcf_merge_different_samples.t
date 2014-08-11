#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 16;
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
push(@output_files, file(@output_subdirs, '1_vcf_merge_different_samples', 'merged.vcf.gz'), file(@output_subdirs, '1_vcf_merge_different_samples', 'merged.vcf.gz.csi'));

ok handle_pipeline(@output_files), 'single step vcf_merge_different_samples pipeline ran and created all expected output files';

my $vcf_file = shift @output_files;
my $samples  = vcf_samples($vcf_file);
is $samples, "SAMPLE01\tSAMPLE02", "merged VCF contains both input samples";
my $scs = VRPipe::StepCmdSummary->get(id => 1)->summary;
my $cmd = VRPipe::Job->get(id => 1)->cmd;
is $scs, 'bcftools merge -O z \@input_vcfs > $output_path && bcftools index $output_path', 'cmd summary shows listed input vcfs';
like $cmd, qr/merge -O z \S+\.vcf.gz \S+\.vcf\.gz /, 'cmd line actually ran with listed vcfs';

# also test the similar vcf_merge_different_samples_control_aware step
($output_dir, $pipeline, $step) = create_single_step_pipeline('vcf_merge_different_samples_control_aware', 'vcf_files');
is_deeply [$step->id, $step->description], [2, 'Merges compressed VCFs using bcftools merge which contain different samples to produce a single VCF containing all the input samples, with the first sample being the one from the VCF tagged with the right control metadata, and identifying that sample as a control in metadata on the output file'], 'vcf_merge_different_samples_control_aware step created and has correct description';

VRPipe::PipelineSetup->create(
    name       => 'vcf_merge_different_samples_control_aware step test',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => file(qw(t data vcf_merge_different_samples_datasource.group.fofn))->absolute,
        options => { metadata_keys => 'group' }
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);

@output_files = ();
for my $de (2 .. 3) {
    @output_subdirs = output_subdirs($de, 2);
    push(@output_files, file(@output_subdirs, '1_vcf_merge_different_samples_control_aware', 'merged.vcf.gz'), file(@output_subdirs, '1_vcf_merge_different_samples_control_aware', 'merged.vcf.gz.csi'));
}

ok handle_pipeline(@output_files), 'single step vcf_merge_different_samples_control_aware pipeline ran and created all expected output files';

my @merged_vcf_files = sort { VRPipe::File->get(path => $a)->meta_value('group') cmp VRPipe::File->get(path => $b)->meta_value('group') } grep { VRPipe::File->get(path => $_)->meta_value('group') } @output_files;
$vcf_file = shift @merged_vcf_files;
$samples  = vcf_samples($vcf_file);
is $samples, "SAMPLE02\tSAMPLE01", "merged VCF contains both input samples in the correct order";
my $control_sample = VRPipe::File->get(path => $vcf_file)->metadata->{sample_control};
is $control_sample, 'foo_SAMPLE02', 'merged VCF has the correct control metadata';
$vcf_file = shift @merged_vcf_files;
$samples  = vcf_samples($vcf_file);
is $samples, "SAMPLE03", "copied VCF contains the input sample";
$control_sample = VRPipe::File->get(path => $vcf_file)->metadata->{sample_control};
is $control_sample, 'bar_SAMPLE03', 'copied VCF has the correct control metadata';

# and test the similar vcf_merge_different_samples_to_indexed_bcf, and that it
# uses the -l option when given lots of input vcfs
($output_dir, $pipeline, $step) = create_single_step_pipeline('vcf_merge_different_samples_to_indexed_bcf', 'vcf_files');
is_deeply [$step->id, $step->description], [3, 'Merges compressed VCFs using bcftools merge which contain different samples to produce a single indexed BCF containing data for all the input samples'], 'vcf_merge_different_samples_to_indexed_bcf step created and has correct description';

VRPipe::PipelineSetup->create(
    name        => 'vcf_merge_different_samples_to_indexed_bcf step test',
    datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'group_all', source => file(qw(t data vcf_merge_different_samples_datasource.30.fofn))->absolute),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);

@output_files = ();
@output_subdirs = output_subdirs(4, 3);
push(@output_files, file(@output_subdirs, '1_vcf_merge_different_samples_to_indexed_bcf', 'merged.bcf'));
push(@output_files, file(@output_subdirs, '1_vcf_merge_different_samples_to_indexed_bcf', 'merged.bcf.csi'));

ok handle_pipeline(@output_files), 'single step vcf_merge_different_samples_to_indexed_bcf pipeline ran and created all expected output files';
$scs = VRPipe::StepCmdSummary->get(id => 2)->summary;
$cmd = VRPipe::Job->get(id => 4)->cmd;
is $scs, 'bcftools merge -O b -l $input_vcfs_fofn > $output_path && bcftools index $output_path', 'cmd summary shows -l input when lots of inputs';
like $cmd, qr/merge -O b -l \S+\.l /, 'cmd line actuall ran with -l option';

finish;

sub vcf_samples {
    my $vcf_file = shift;
    open(my $fh, "zcat $vcf_file |");
    my $samples;
    while (<$fh>) {
        next unless /^#CHROM.+FORMAT\s+(.+)$/;
        $samples = $1;
        last;
    }
    close($fh);
    return $samples;
}
