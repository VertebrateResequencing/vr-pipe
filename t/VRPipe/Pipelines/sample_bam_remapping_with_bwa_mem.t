#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => 'VRPIPE_TEST_PIPELINES',
        required_exe => [qw(bwa samtools bamtofastq bamstreamingmarkduplicates bamsort)]
    );
    use TestPipelines;
}

my $mapping_output_dir = get_output_dir('sample_bam_remapping_with_bwa_mem');

ok my $mapping_pipeline = VRPipe::Pipeline->create(name => 'sample_bam_remapping_with_bwa_mem'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($mapping_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(fasta_index sequence_dictionary bwa_index samtools_split_by_readgroup bam_metadata bamtofastq bwa_mem_to_bam biobambam_bammerge_and_streaming_mark_duplicates)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($mapping_output_dir, 'ref');
$mapping_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $mapping_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'sample_bam_remapping_with_bwa_mem',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => file(qw(t data sample_remapping.fofnwm))->absolute->stringify,
        options => { metadata_keys => 'sample' }
    ),
    output_root => $mapping_output_dir,
    pipeline    => $mapping_pipeline,
    options     => {
        reference_fasta            => $ref_fa,
        reference_assembly_name    => 'SSuis1',
        reference_public_url       => 'ftp://s.suis.com/ref.fa',
        reference_species          => 'S.Suis',
        sequence_dictionary_memory => 150,
        sequence_dictionary_time   => 1,
        bwa_index_options          => '-a is',
        bwa_index_memory           => 150,
        bwa_index_time             => 1,
        store_original_pg_chain    => 0,
        fastq_chunk_size           => 8000,
        rescue_orphans             => 1,
        bwa_mem_post_processing    => 'bamsort inputformat=sam fixmates=1 adddupmarksupport=1',
        samtools_merge_options     => '-cpu',
        cleanup                    => 0
    }
);
ok handle_pipeline(), 'pipeline ran ok';

finish;
