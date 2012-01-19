#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],
                    required_exe => [qw(bamcheck)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::bam_metadata');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_metadata', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Takes a bam file and associates metadata with the file in the VRPipe database, making the bam file usable in other bam-related Steps'], 'bam_metadata step created and has correct description';

# test as part of a pipeline
my $setup = VRPipe::PipelineSetup->get(name => 'bm_setup',
                                       datasource => VRPipe::DataSource->get(type => 'fofn', method => 'all', source => file(qw(t data datasource.bam_fofn))->absolute),
                                       output_root => $output_dir,
                                       pipeline => $pipeline,
                                       options => {});
ok handle_pipeline(), 'pipeline ran ok with no options';

my $bam = VRPipe::File->get(path => file(qw(t data remapping_bams 8324_8.pe.bam))->absolute);
my $meta = $bam->metadata;
is_deeply [$meta->{reads}, $meta->{paired}, $meta->{forward_reads}, $meta->{reads_mapped}, $meta->{sd_insert_size}, $meta->{rmdup_reads}], [500, 1, 250, 233, 66.8, undef], 'metadata was correctly applied to a bam file';

$setup = VRPipe::PipelineSetup->get(name => 'bm_setup_with_d',
                                    datasource => VRPipe::DataSource->get(type => 'fofn', method => 'all', source => file(qw(t data datasource.bam_fofn))->absolute),
                                    output_root => $output_dir,
                                    pipeline => $pipeline,
                                    options => {bamcheck_options => '-d'});
ok handle_pipeline(), 'pipeline ran ok with bamcheck_options';

is $bam->metadata->{rmdup_reads}, 500, 'rerunning with bamcheck_options => -d added the rmdup_reads metadata';

finish;
