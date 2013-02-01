#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES GATK2)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::gatk_reduce_reads');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('gatk_reduce_reads', 'bam_files');
VRPipe::StepAdaptor->create(pipeline => $pipeline, to_step => 1, adaptor_hash => { 'bai_files' => { data_element => 0 }, 'bam_files' => { data_element => 0 } });
is_deeply [$step->id, $step->description], [1, 'Applies GATK ReduceReads to one or more BAM files which reduces the BAM file using read based compression that keeps only essential information for variant calling'], 'gatk_reduce_reads step created and has correct description';

my $ref_fa = file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify;

# test as part of a pipeline
my $setup = VRPipe::PipelineSetup->create(
    name        => 'reduce_reads_setup',
    datasource  => VRPipe::DataSource->create(type => 'fofn_with_metadata', method => 'grouped_by_metadata', source => file(qw(t data hs_chr20.qc.bambai.fofnwmd))->absolute, options => { metadata_keys => 'lane' }),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options => { reference_fasta => $ref_fa }
);

# order of elements not stable with grouped_by_metadata? figure out the order
# now
$setup->datasource->elements;
my %bams;
foreach my $e (VRPipe::DataElement->search({})) {
    my $lane = $e->result->{group};
    $bams{"hs_chr20.$lane.bam"} = $e->id;
}

my @output_files;
while (my ($bam, $element_id) = each %bams) {
    my @output_subdirs = output_subdirs($element_id);
    $bam =~ s/bam$/reduced.bam/;
    push(@output_files, file(@output_subdirs, '1_gatk_reduce_reads', $bam));
}

ok handle_pipeline(@output_files), 'gatk_reduce_reads single-step pipeline ran ok and created expected output files';

finish;
