#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(samtools smalt)]
    );
    
    use TestPipelines;
}

my $mapping_output_dir = get_output_dir('smalt_mapping');

ok my $mapping_pipeline = VRPipe::Pipeline->create(name => 'fastq_mapping_with_smalt'), 'able to get a pre-written pipeline';
my @s_names;
foreach my $stepmember ($mapping_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(sequence_dictionary smalt_index fastq_import fastq_metadata fastq_split fastq_decompress smalt_map_to_sam sam_to_fixed_bam bam_add_readgroup bam_merge_lane_splits)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($mapping_output_dir, 'ref');
$mapping_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);
my $mapping_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 's_suis mapping',
    datasource => VRPipe::DataSource->create(
        type    => 'sequence_index',
        method  => 'lane_fastqs',
        source  => file(qw(t data datasource.sequence_index)),
        options => { local_root_dir => dir(".")->absolute->stringify }
    ),
    output_root => $mapping_output_dir,
    pipeline    => $mapping_pipeline,
    options     => {
        fastq_chunk_size             => 8000,
        reference_fasta              => $ref_fa,
        reference_assembly_name      => 'SSuis1',
        reference_public_url         => 'ftp://s.suis.com/ref.fa',
        reference_species            => 'S.Suis',
        smalt_index_options          => '-k 13 -s 4',
        fixed_bam_seq_from_reference => 1,
        cleanup                      => 0,
        sequence_dictionary_memory   => 150,
        sequence_dictionary_time     => 1,
        smalt_index_memory           => 150,
        smalt_index_time             => 1,
        fastq_import_memory          => 150,
        fastq_import_time            => 1,
        fastq_metadata_memory        => 150,
        fastq_metadata_time          => 1,
        fastq_split_memory           => 150,
        fastq_split_time             => 1,
        sam_to_fixed_bam_memory      => 150,
        sam_to_fixed_bam_time        => 1,
        bam_merge_lane_splits_memory => 150,
        bam_merge_lane_splits_time   => 1
    }
);

ok handle_pipeline(), 'pipeline ran ok';

is_deeply [VRPipe::StepState->create(pipelinesetup => 1, stepmember => 2, dataelement => 1)->cmd_summary->summary, VRPipe::StepState->create(pipelinesetup => 1, stepmember => 7, dataelement => 1)->cmd_summary->summary, VRPipe::StepState->create(pipelinesetup => 1, stepmember => 8, dataelement => 1)->cmd_summary->summary, VRPipe::StepState->create(pipelinesetup => 1, stepmember => 9, dataelement => 1)->cmd_summary->summary], ['smalt index -k 13 -s 4 $index_base $reference_fasta', 'smalt map -f samsoft -i 400 -o $sam_file $index_base $fastq_file(s)', 'samtools view -bSu $sam_file -T $reference_fasta | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - $reference_fasta > $fixed_bam_file', 'java $jvm_args -jar AddOrReplaceReadGroups.jar INPUT=$bam_file OUTPUT=$rg_added_bam_file RGID=$lane RGLB=$library RGPL=$platform RGPU=$platform_unit RGSM=$sample RGCN=$centre RGDS=$study VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0'], 'cmd summaries for the major steps were as expected';

finish;
