#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(samtools)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::bam_processing');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_processing', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Runs the input BAM or CRAM through command line(s) given by the user'], 'bam_processing step created and has correct description';

# test as part of a pipeline
my $ds = VRPipe::DataSource->create(type => 'fofn_with_metadata', method => 'all', source => file(qw(t data sample_remapping.fofnwm))->absolute);
my $setup = VRPipe::PipelineSetup->create(
    name        => 'processing_setup',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        command_line => '$samtools view -l ${library} -u -o $output_bam $input_bam',
        index_output => 1,
    }
);
ok handle_pipeline(), 'pipeline ran ok with no options';

finish;
