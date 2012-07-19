#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::chunk_genomic_region');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('chunk_genomic_region', 'no_key');
is_deeply [$step->id, $step->description], [1, 'Generate a chromosomal regions file, split according to chunk size, from fasta reference index file or specific regions file'], 'chunk_genomic_region step created and has correct description';

my $genomic_region_file  = VRPipe::File->create(path => file(qw(t data human_g1k_v37.fasta.fai))->absolute)->path;
my $ploidy_file          = VRPipe::File->create(path => file(qw(t data ploidy_definition))->absolute)->path;
my $chunked_regions_file = VRPipe::File->create(path => file($output_dir, 'chunked_regions_file.txt'))->path;
my $ok = VRPipe::Steps::chunk_genomic_region->write_chunked_regions_file($genomic_region_file->stringify, $chunked_regions_file->stringify, '1 2 11 X Y MT', 1000000, 0, $ploidy_file->stringify);

my $ploidy_def = do $ploidy_file;
my $ploidy = { default => 2,
               X       => [{ from => 1, to => 60_000, M => 1 }, { from => 2_699_521, to => 154_931_043, M => 1 },],
               Y  => [{ from => 1, to => 59_373_566, M => 1, F => 0 },],
               MT => [{ from => 1, to => 16_569,     M => 1, F => 1 },], };
is_deeply [$ploidy_def], [$ploidy], 'ploidy file read correctly';

# test as part of a pipeline
my $setup = VRPipe::PipelineSetup->create(name        => 'chunk_genomic_region_setup',
                                          datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data improvement_datasource.fofn))->absolute),
                                          output_root => $output_dir,
                                          pipeline    => $pipeline,
                                          options     => {
                                                       genomic_region_file => $genomic_region_file->stringify,
                                                       chrom_list          => '1 2 11 20 X Y MT',
                                                       chunk_size          => '10000000',
                                                       ploidy              => $ploidy_file->stringify });

ok handle_pipeline(), 'single-step pipeline ran ok';

finish;
