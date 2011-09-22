#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPipelines;
}

my $mapping_output_dir = get_output_dir('mapping_with_improvement');

ok my $mapping_pipeline = VRPipe::Pipeline->get(name => '1000genomes_illumina_mapping_with_improvement'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($mapping_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(sequence_dictionary
                         bwa_index
                         fastq_import
                         fastq_metadata
                         fastq_split
                         bwa_aln_fastq
                         bwa_sam
                         sam_to_fixed_bam
                         bam_merge_lane_splits
                         bam_realignment_around_known_indels
                         bam_fix_mates
                         bam_recalibrate_quality_scores
                         bam_calculate_bq
                         bam_reheader
                         bam_stats)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($mapping_output_dir, 'ref');
$mapping_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);
my $mapping_pipelinesetup = VRPipe::PipelineSetup->get(name => 's_suis mapping',
                                                       datasource => VRPipe::DataSource->get(type => 'sequence_index',
                                                                                             method => 'lane_fastqs',
                                                                                             source => file(qw(t data datasource.sequence_index)),
                                                                                             options => { local_root_dir => dir(".")->absolute->stringify }),
                                                       output_root => $mapping_output_dir,
                                                       pipeline => $mapping_pipeline,
                                                       options => {fastq_chunk_size => 8000,
                                                                   reference_fasta => $ref_fa,
                                                                   reference_assembly_name => 'SSuis1',
                                                                   reference_public_url => 'ftp://s.suis.com/ref.fa',
                                                                   reference_species => 'S.Suis',
                                                                   bwa_index_options => '-a is',
                                                                   cleanup => 0,
                                                                   sequence_dictionary_memory => 150,
                                                                   sequence_dictionary_time => 1,
                                                                   bwa_index_memory => 150,
                                                                   bwa_index_time => 1,
                                                                   fastq_import_memory => 150,
                                                                   fastq_import_time => 1,
                                                                   fastq_metadata_memory => 150,
                                                                   fastq_metadata_time => 1,
                                                                   fastq_split_memory => 150,
                                                                   fastq_split_time => 1,
                                                                   bwa_aln_fastq_memory => 150,
                                                                   bwa_aln_fastq_time => 1,
                                                                   bwa_sam_memory => 150,
                                                                   bwa_sam_time => 1,
                                                                   sam_to_fixed_bam_memory => 150,
                                                                   sam_to_fixed_bam_time => 1,
                                                                   bam_merge_lane_splits_memory => 150,
                                                                   bam_merge_lane_splits_time => 1});

handle_pipeline();

is_deeply [VRPipe::StepState->get(pipelinesetup => 1, stepmember => 2, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 6, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 7, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 8, dataelement => 1)->cmd_summary->summary],
          ['bwa index -a is $reference_fasta',
           'bwa aln -q 15 -f $sai_file $reference_fasta $fastq_file',
           'bwa sampe -a 600 -r $rg_line -f $sam_file $reference_fasta $sai_file(s) $fastq_file(s)',
           'samtools view -bSu $sam_file | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - $reference_fasta > $fixed_bam_file'],
          'cmd summaries for the major steps were as expected';

finish;