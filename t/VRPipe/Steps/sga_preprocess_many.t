#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(sga)]
    );
    use TestPipelines;
    use_ok('VRPipe::Steps::sga_preprocess_many');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('sga_preprocess_many', 'fastq_files');
is_deeply [$step->id, $step->description], [1, 'Prepare fastq files for assembly with sga.'], 'sga_preprocess_many step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'sga_preprocess_many',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => file(qw(t data fq.fofn_with_metadata))->absolute->stringify,
        options => {
            metadata_keys => 'lane',
        }
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        sga_preprocess_options => '--min-length=50',
    }
);

ok handle_pipeline(), 'sga_preprocess_many pipeline ran ok';

# my @output_subdirs = output_subdirs(1);
# ok handle_pipeline(file(@output_subdirs, '1_sga_preprocess_many', "dbsnp_132.b37.chr20_reduced.aln.0.vcf.gz")), 'gatk_left_align_variants pipeline ran ok, generating the expected file';
