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

my $mapping_output_dir = get_output_dir('fastq_mapping_with_bwa_mem');

ok my $mapping_pipeline = VRPipe::Pipeline->create(name => 'fastq_mapping_with_bwa_mem'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($mapping_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(fasta_index sequence_dictionary bwa_index fastq_import fastq_metadata fastq_split bwa_mem_fastq sam_to_fixed_bam bam_merge_lane_splits bam_index)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($mapping_output_dir, 'ref');
$mapping_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);
my $mapping_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'fastq_mapping_with_bwa_mem',
    datasource => VRPipe::DataSource->create(
        type    => 'sequence_index',
        method  => 'lane_fastqs',
        source  => file(qw(t data datasource.sequence_index)),
        options => { local_root_dir => dir(".")->absolute->stringify }
    ),
    output_root => $mapping_output_dir,
    pipeline    => $mapping_pipeline,
    options     => {
        reference_fasta               => $ref_fa,
        reference_assembly_name       => 'SSuis1',
        reference_public_url          => 'ftp://s.suis.com/ref.fa',
        reference_species             => 'S.Suis',
        bwa_index_options             => '-a is',
        fastq_chunk_size              => 8000,
        uncompressed_fixed_bam_output => 0,
        sequence_dictionary_memory    => 150,
        sequence_dictionary_time      => 1,
        bwa_index_memory              => 150,
        bwa_index_time                => 1,
        bam_metadata_memory           => 150,
        bam_metadata_time             => 1,
        fastq_split_memory            => 150,
        fastq_split_time              => 1,
        bwa_mem_fastq_memory          => 150,
        bwa_mem_fastq_time            => 1,
        sam_to_fixed_bam_memory       => 150,
        sam_to_fixed_bam_time         => 1,
        bam_merge_lane_splits_memory  => 150,
        bam_merge_lane_splits_time    => 1,
        cleanup                       => 1
    }
);
ok handle_pipeline(), 'pipeline ran ok';

#*** needs proper tests...

finish;
