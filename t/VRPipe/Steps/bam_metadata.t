#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(bamcheck)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::bam_metadata');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_metadata', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Takes a bam file and associates metadata with the file in the VRPipe database, making the bam file usable in other bam-related Steps'], 'bam_metadata step created and has correct description';

# test as part of a pipeline
my $ds = VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data datasource.bam_fofn))->absolute);
my $setup = VRPipe::PipelineSetup->create(
    name        => 'bm_setup',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);
ok handle_pipeline(), 'pipeline ran ok with no options';

my $bam = VRPipe::File->create(path => file(qw(t data remapping_bams 8324_8.pe.bam))->absolute);
my $meta = $bam->metadata;
is_deeply [$meta->{reads}, $meta->{paired}, $meta->{forward_reads}, $meta->{reads_mapped}, $meta->{sd_insert_size}, $meta->{rmdup_reads}], [500, 1, 250, 233, 66.8, 500], 'metadata was correctly applied to a bam file';

($output_dir, $pipeline, $step) = create_single_step_pipeline('bamcheck', 'bam_files');
$setup = VRPipe::PipelineSetup->create(
    name        => 'bm_setup_with_d',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta  => file(qw(t data S_suis_P17.fa))->absolute->stringify,
        bamcheck_options => '-d'
    }
);
ok handle_pipeline(), 'bamcheck pipeline ran ok with bamcheck_options';

$bam->reselect_values_from_db;
is $bam->metadata->{rmdup_bases_mapped_c}, 12289, 'running bamcheck with bamcheck_options => -d added the rmdup_bases_mapped_c metadata';

finish;
