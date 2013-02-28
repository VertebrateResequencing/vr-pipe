#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES )],
        required_exe => [qw(samtools)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('run_bam2fastq');

ok my $pipeline = VRPipe::Pipeline->create(name => 'run_bam2fastq'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(run_bam2fastq)], 'the pipeline has the correct steps';

# pipeline requires certain metadata on the input bams, so we just manually set that now
VRPipe::File->create(
    path     => file(qw(t data 2822_6.pe.bam))->absolute,
    metadata => {
        lane          => '2822_6',
        reads         => 400,
        bases         => 23000,
        forward_reads => 200,
        reverse_reads => 200,
        paired        => 1
    }
);
VRPipe::File->create(
    path     => file(qw(t data 2822_6.improved.pe.bam))->absolute,
    metadata => {
        lane          => '2822_6',
        reads         => 400,
        bases         => 23000,
        forward_reads => 200,
        reverse_reads => 200,
        paired        => 1
    }
);

my $pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'run_bam2fastq',
    datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data datasource.btf_fofn))->absolute),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        bam2fastq_exe => '/lustre/scratch106/user/cj5/bam2fastq-jts/bam2fastq',
    }
);

ok handle_pipeline(), 'pipeline ran and created all expected output files';

finish;
