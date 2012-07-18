#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],
                    required_exe => [qw(bamcheck)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::bam_metadata');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_metadata_with_sex', 'bam_files');

my $ds = VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data datasource.bam_fofn))->absolute);
my $setup = VRPipe::PipelineSetup->create(name        => 'bm_setup',
                                          datasource  => $ds,
                                          output_root => $output_dir,
                                          pipeline    => $pipeline,
                                          options     => { sample_sex_file => file(qw(t data datasource.sample.sex))->absolute->stringify, });
ok handle_pipeline(), 'pipeline ran ok';

my $bam = VRPipe::File->get(path => file(qw(t data remapping_bams 8324_8.pe.bam))->absolute);
my $meta = $bam->metadata;
is_deeply [$meta->{reads}, $meta->{paired}, $meta->{forward_reads}, $meta->{reads_mapped}, $meta->{sd_insert_size}, $meta->{rmdup_reads}, $meta->{sex}], [500, 1, 250, 233, 66.8, 500, 'M'], 'metadata including sex was correctly applied to a bam file';

finish;
