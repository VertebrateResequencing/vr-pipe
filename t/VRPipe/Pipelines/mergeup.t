#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPipelines;
}

ok my $mapping_pipeline = VRPipe::Pipeline->get(name => 'fastq_mapping_with_bwa'), 'able to get the fastq_mapping_with_bwa pipeline';
ok my $merge_lanes_pipeline = VRPipe::Pipeline->get(name => 'merge_lanes'), 'able to get the merge_lanes pipeline';
ok my $merge_libraries_pipeline = VRPipe::Pipeline->get(name => 'merge_libraries_and_split'), 'able to get the merge_libraries_and_split pipeline';

my $mapping_dir = get_output_dir('mapping_mergeup_test');

my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($mapping_dir, 'ref');
$mapping_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $mapping_pipelinesetup = VRPipe::PipelineSetup->get(name => 's_suis mapping',
                                                       datasource => VRPipe::DataSource->get(type => 'sequence_index',
                                                                                             method => 'lane_fastqs',
                                                                                             source => file(qw(t data datasource.sequence_index)),
                                                                                             options => { local_root_dir => dir(".")->absolute->stringify }),
                                                       output_root => $mapping_dir,
                                                       pipeline => $mapping_pipeline,
                                                       options => {fastq_chunk_size => 8000,
                                                                   reference_fasta => $ref_fa,
                                                                   reference_assembly_name => 'SSuis1',
                                                                   reference_public_url => 'ftp://s.suis.com/ref.fa',
                                                                   reference_species => 'S.Suis',
                                                                   bwa_index_options => '-a is',
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

my $build_dir = get_output_dir('build_test');

my $merge_lanes_pipelinesetup = VRPipe::PipelineSetup->get(name => 's_suis merge lanes',
                                                           datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                                                 method => 'group_by_metadata',
                                                                                                 source => 's_suis mapping[9:merged_lane_bams]',
                                                                                                 options => { metadata_keys => 'sample|platform|library' } ),
                                                           output_root => $build_dir,
                                                           pipeline => $merge_lanes_pipeline,
                                                           options => { bam_tags_to_strip => 'OQ|XM|XG|XO',
                                                                        bam_merge_keep_single_paired_separate => 1,
                                                                        bam_merge_memory => 200,
                                                                        cleanup => 1 });

my $merge_libraries_pipelinesetup = VRPipe::PipelineSetup->get(name => 's_suis merge libraries',
                                                               datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                                                     method => 'group_by_metadata',
                                                                                                     source => 's_suis merge lanes[3:markdup_bam_files]',
                                                                                                     options => { metadata_keys => 'sample|platform' } ),
                                                               output_root => $build_dir,
                                                               pipeline => $merge_libraries_pipeline,
                                                               options => { bam_merge_keep_single_paired_separate => 0,
                                                                            reference_fasta => $ref_fa,
                                                                            reference_assembly_name => 'SSuis1',
                                                                            reference_public_url => 'ftp://s.suis.com/ref.fa',
                                                                            reference_species => 'S.Suis',
                                                                            bam_merge_memory => 200,
                                                                            cleanup => 1 });

ok handle_pipeline(), 'pipelines ran ok';

is_deeply [VRPipe::StepState->get(pipelinesetup => 1, stepmember => 2, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 6, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 7, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 8, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 2, stepmember => 12, dataelement => 5)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 2, stepmember => 13, dataelement => 5)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 3, stepmember => 14, dataelement => 8)->cmd_summary->summary],
          ['bwa index -a is $reference_fasta',
           'bwa aln -q 15 -f $sai_file $reference_fasta $fastq_file',
           'bwa sampe -a 600 -r $rg_line -f $sam_file $reference_fasta $sai_file(s) $fastq_file(s)',
           'samtools view -bSu $sam_file | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - $reference_fasta > $fixed_bam_file',
           'java $jvm_args -jar MergeSamFiles.jar INPUT=$bam_file(s) OUTPUT=$merged_bam VALIDATION_STRINGENCY=SILENT',
           'java $jvm_args -jar MarkDuplicates.jar INPUT=$bam_file OUTPUT=$markdup_bam_file ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT',
           'java $jvm_args -jar MergeSamFiles.jar INPUT=$bam_file(s) OUTPUT=$merged_bam VALIDATION_STRINGENCY=SILENT'],
          'cmd summaries for the major steps were as expected';

# my @expected_header_lines = ("\@HD\tVN:1.0\tSO:coordinate",
#                              "\@SQ\tSN:fake_chr1\tLN:290640\tM5:55f9584cf1f4194f13bbdc0167e0a05f\tUR:ftp://s.suis.com/ref.fa\tAS:SSuis1\tSP:S.Suis",
#                              "\@SQ\tSN:fake_chr2\tLN:1716851\tM5:6dd2836053e5c4bd14ad49b5b2f2eb88\tUR:ftp://s.suis.com/ref.fa\tAS:SSuis1\tSP:S.Suis",
#                              "\@RG\tID:2822_6\tLB:LIB01\tSM:SAMPLE01\tPI:200\tCN:SC\tPL:ILLUMINA\tDS:STUDY01",
#                              "\@PG\tID:bwa_index\tPN:bwa\tVN:0.5.9-r16\tCL:bwa index -a is \$reference_fasta",
#                              "\@PG\tID:bwa_aln_fastq\tPN:bwa\tPP:bwa_index\tVN:0.5.9-r16\tCL:bwa aln -q 15 -f \$sai_file \$reference_fasta \$fastq_file",
#                              "\@PG\tID:bwa_sam\tPN:bwa\tPP:bwa_aln_fastq\tVN:0.5.9-r16\tCL:bwa sampe -a 600 -r \$rg_line -f \$sam_file \$reference_fasta \$sai_file(s) \$fastq_file(s)",
#                              "\@PG\tID:sam_to_fixed_bam\tPN:samtools\tPP:bwa_sam\tVN:0.1.17 (r973:277)\tCL:samtools view -bSu \$sam_file | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - \$reference_fasta > \$fixed_bam_file");
# 
# my @header_lines = get_bam_header($test_bam);
# is_deeply \@header_lines, \@expected_header_lines, 'bam header is okay';

finish;