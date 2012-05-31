#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 18;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES PICARD GATK)],
                    required_exe => [qw(samtools bwa)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::picard');
}

my $picard = VRPipe::Steps::picard->new(picard_path => $ENV{PICARD});
my $picard_version = $picard->determine_picard_version();
my $bwa_version = VRPipe::StepCmdSummary->determine_version('bwa', '^Version: (.+)$');
my $samtools_version = VRPipe::StepCmdSummary->determine_version('samtools', '^Version: (.+)$');

ok my $mapping_pipeline = VRPipe::Pipeline->get(name => 'fastq_mapping_with_bwa'), 'able to get the fastq_mapping_with_bwa pipeline';
ok my $merge_lanes_pipeline = VRPipe::Pipeline->get(name => 'bam_merge_lanes'), 'able to get the bam_merge_lanes pipeline';
ok my $merge_libraries_pipeline = VRPipe::Pipeline->get(name => 'bam_merge_1000_genomes_libraries_with_split'), 'able to get the bam_merge_1000_genomes_libraries_with_split pipeline';
ok my $release_pipeline = VRPipe::Pipeline->get(name => '1000genomes_release'), 'able to get the 1000genomes_release pipeline';

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
                                                                                                 options => { metadata_keys => 'analysis_group|population|sample|platform|library' } ),
                                                           output_root => $build_dir,
                                                           pipeline => $merge_lanes_pipeline,
                                                           options => { bam_tags_to_strip => 'OQ XM XG XO',
                                                                        bam_merge_keep_single_paired_separate => 1,
                                                                        bam_merge_memory => 200,
                                                                        cleanup => 1 });

my $merge_libraries_pipelinesetup = VRPipe::PipelineSetup->get(name => 's_suis merge libraries',
                                                               datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                                                     method => 'group_by_metadata',
                                                                                                     source => 's_suis merge lanes[3:markdup_bam_files]',
                                                                                                     options => { metadata_keys => 'analysis_group|population|sample|platform' } ),
                                                               output_root => $build_dir,
                                                               pipeline => $merge_libraries_pipeline,
                                                               options => { bam_merge_keep_single_paired_separate => 0,
                                                                            reference_fasta => $ref_fa,
                                                                            reference_assembly_name => 'SSuis1',
                                                                            reference_public_url => 'ftp://s.suis.com/ref.fa',
                                                                            reference_species => 'S.Suis',
                                                                            bam_merge_memory => 200,
                                                                            split_bam_make_unmapped => 1,
                                                                            cleanup => 1,
                                                                            delete_input_bams => 1,
                                                                            remove_merged_bams => 1 });

my $release_pipeline_setup = VRPipe::PipelineSetup->get(name => 's_suis release',
                                                               datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                                                     method => 'all',
                                                                                                     source => 's_suis merge libraries[4:split_bam_files]',
                                                                                                     options => { filter => 'split_sequence#^(fake_chr2|unmapped)$',
                                                                                                                  filter_after_grouping => 0 } ),
                                                               output_root => $build_dir,
                                                               pipeline => $release_pipeline,
                                                               options => { release_date => '19790320',
                                                                            sequence_index => file(qw(t data datasource.sequence_index))->absolute->stringify,
                                                                            rg_from_pu => 0 });

my @mapping_files;
my %bams = ('2822_6.pe.bam' => 1, '2822_6.se.bam' => 1, '2822_7.pe.bam' => 2, '2823_4.pe.bam' => 3, '8324_8.pe.bam' => 4);
while (my ($bam, $element_id) = each %bams) {
    my @output_subdirs = output_subdirs($element_id);
    push(@mapping_files, file(@output_subdirs, '9_bam_merge_lane_splits', $bam));
}


handle_pipeline(@mapping_files);
# linked pipelines only create their dataelements after the previous pipeline
# has finished, so only now can we figure out remaining expected output file
# paths

my @split_files;
my @split_files_removed;
foreach my $element (@{$merge_libraries_pipelinesetup->datasource->elements}) {
    my @output_subdirs = output_subdirs($element->id, 3);
    foreach my $file ('fake_chr1.pe.1.bam', 'fake_chr2.pe.1.bam', 'unmapped.pe.1.bam') {
        push(@split_files, file(@output_subdirs, '4_bam_split_by_sequence', $file));
    }
    if ($element->result->{group} =~ /SAMPLE01/) {
        push @split_files_removed, (1,1,1);
    } else {
        push @split_files_removed, (0,0,0);
    }
}

my @release_files;
my @release_files_removed;
foreach my $element (@{$release_pipeline_setup->datasource->elements}) {
    my @output_subdirs = output_subdirs($element->id, 4);
    foreach my $file ('fake_chr2.pe.1.bam', 'unmapped.pe.1.bam') {
        push(@release_files, file(@output_subdirs, '1_dcc_metadata', $file));
        push(@release_files, file(@output_subdirs, '1_dcc_metadata', $file.'.bai'));
        push(@release_files, file(@output_subdirs, '3_bam_stats', $file.'.bas'));
        push(@release_files, file(@output_subdirs, '4_md5_file_production', $file.'.md5'));
        push(@release_files, file(@output_subdirs, '5_md5_file_production', $file.'.bai.md5'));
        push(@release_files, file(@output_subdirs, '6_md5_file_production', $file.'.bas.md5'));
    }
    if ($element->result->{group} =~ /SAMPLE01/) {
        push @release_files_removed, (1,1,1,1,1,1,1,1,1,1,1,1);
    } else {
        push @release_files_removed, (0,0,0,0,0,0,0,0,0,0,0,0);
    }
}

ok handle_pipeline(@mapping_files, @split_files, @release_files), 'pipelines ran ok and correct output files created';

is_deeply [VRPipe::StepState->get(pipelinesetup => 1, stepmember => 2, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 6, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 7, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 8, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 2, stepmember => 11, dataelement => 5)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 2, stepmember => 12, dataelement => 5)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 3, stepmember => 13, dataelement => 8)->cmd_summary->summary],
          ['bwa index -a is $reference_fasta',
           'bwa aln -q 15 -f $sai_file $reference_fasta $fastq_file',
           'bwa sampe -a 600 -r $rg_line -f $sam_file $reference_fasta $sai_file(s) $fastq_file(s)',
           'samtools view -bSu $sam_file | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - $reference_fasta > $fixed_bam_file',
           'java $jvm_args -jar MergeSamFiles.jar INPUT=$bam_file(s) OUTPUT=$merged_bam VALIDATION_STRINGENCY=SILENT',
           'java $jvm_args -jar MarkDuplicates.jar INPUT=$bam_file OUTPUT=$markdup_bam_file ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT',
           'java $jvm_args -jar MergeSamFiles.jar INPUT=$bam_file(s) OUTPUT=$merged_bam VALIDATION_STRINGENCY=SILENT'],
          'cmd summaries for the major steps were as expected';

my @expected_header_lines = ("\@HD\tVN:1.0\tSO:coordinate",
                             "\@SQ\tSN:fake_chr1\tLN:290640\tM5:55f9584cf1f4194f13bbdc0167e0a05f\tUR:ftp://s.suis.com/ref.fa\tAS:SSuis1\tSP:S.Suis",
                             "\@SQ\tSN:fake_chr2\tLN:1716851\tM5:6dd2836053e5c4bd14ad49b5b2f2eb88\tUR:ftp://s.suis.com/ref.fa\tAS:SSuis1\tSP:S.Suis",
                             "\@RG\tID:2822_6\tLB:LIB01\tSM:SAMPLE01\tPI:200\tCN:SC\tPL:ILLUMINA\tDS:STUDY01",
                             "\@RG\tID:2822_7\tLB:LIB01\tSM:SAMPLE01\tPI:205\tCN:SC\tPL:ILLUMINA\tDS:STUDY01",
                             "\@RG\tID:2823_4\tLB:LIB02\tSM:SAMPLE01\tPI:195\tCN:SC\tPL:ILLUMINA\tDS:STUDY01",
                             "\@PG\tID:bwa_index\tPN:bwa\tVN:$bwa_version\tCL:bwa index -a is \$reference_fasta",
                             "\@PG\tID:bwa_aln_fastq\tPN:bwa\tPP:bwa_index\tVN:$bwa_version\tCL:bwa aln -q 15 -f \$sai_file \$reference_fasta \$fastq_file",
                             "\@PG\tID:bwa_sam\tPN:bwa\tPP:bwa_aln_fastq\tVN:$bwa_version\tCL:bwa sampe -a 600 -r \$rg_line -f \$sam_file \$reference_fasta \$sai_file(s) \$fastq_file(s)",
                             "\@PG\tID:bwa_sam.1\tPN:bwa\tPP:bwa_aln_fastq\tVN:$bwa_version\tCL:bwa sampe -a 615 -r \$rg_line -f \$sam_file \$reference_fasta \$sai_file(s) \$fastq_file(s)",
                             "\@PG\tID:bwa_sam.2\tPN:bwa\tPP:bwa_aln_fastq\tVN:$bwa_version\tCL:bwa sampe -a 585 -r \$rg_line -f \$sam_file \$reference_fasta \$sai_file(s) \$fastq_file(s)",
                             "\@PG\tID:sam_to_fixed_bam\tPN:samtools\tPP:bwa_sam\tVN:$samtools_version\tCL:samtools view -bSu \$sam_file | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - \$reference_fasta > \$fixed_bam_file",
                             "\@PG\tID:sam_to_fixed_bam.1\tPN:samtools\tPP:bwa_sam.1\tVN:$samtools_version\tCL:samtools view -bSu \$sam_file | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - \$reference_fasta > \$fixed_bam_file",
                             "\@PG\tID:sam_to_fixed_bam.2\tPN:samtools\tPP:bwa_sam.2\tVN:$samtools_version\tCL:samtools view -bSu \$sam_file | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - \$reference_fasta > \$fixed_bam_file",
                             "\@PG\tID:bam_merge\tPN:picard\tPP:sam_to_fixed_bam\tVN:$picard_version\tCL:java \$jvm_args -jar MergeSamFiles.jar INPUT=\$bam_file(s) OUTPUT=\$merged_bam VALIDATION_STRINGENCY=SILENT",
                             "\@PG\tID:bam_merge.1\tPN:picard\tPP:sam_to_fixed_bam.1\tVN:$picard_version\tCL:java \$jvm_args -jar MergeSamFiles.jar INPUT=\$bam_file(s) OUTPUT=\$merged_bam VALIDATION_STRINGENCY=SILENT",
                             "\@PG\tID:bam_merge.2\tPN:picard\tPP:sam_to_fixed_bam.2\tVN:$picard_version\tCL:java \$jvm_args -jar MergeSamFiles.jar INPUT=\$bam_file(s) OUTPUT=\$merged_bam VALIDATION_STRINGENCY=SILENT",
                             "\@PG\tID:bam_mark_duplicates\tPN:picard\tPP:bam_merge\tVN:$picard_version\tCL:java \$jvm_args -jar MarkDuplicates.jar INPUT=\$bam_file OUTPUT=\$markdup_bam_file ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT",
                             "\@PG\tID:bam_mark_duplicates.1\tPN:picard\tPP:bam_merge.1\tVN:$picard_version\tCL:java \$jvm_args -jar MarkDuplicates.jar INPUT=\$bam_file OUTPUT=\$markdup_bam_file ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT",
                             "\@PG\tID:bam_mark_duplicates.2\tPN:picard\tPP:bam_merge.2\tVN:$picard_version\tCL:java \$jvm_args -jar MarkDuplicates.jar INPUT=\$bam_file OUTPUT=\$markdup_bam_file ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT",
                             "\@PG\tID:bam_merge.3\tPN:picard\tPP:bam_mark_duplicates\tVN:$picard_version\tCL:java \$jvm_args -jar MergeSamFiles.jar INPUT=\$bam_file(s) OUTPUT=\$merged_bam VALIDATION_STRINGENCY=SILENT",
                             "\@PG\tID:bam_merge.1.3\tPN:picard\tPP:bam_mark_duplicates.1\tVN:$picard_version\tCL:java \$jvm_args -jar MergeSamFiles.jar INPUT=\$bam_file(s) OUTPUT=\$merged_bam VALIDATION_STRINGENCY=SILENT",
                             "\@PG\tID:bam_merge.2.3\tPN:picard\tPP:bam_mark_duplicates.2\tVN:$picard_version\tCL:java \$jvm_args -jar MergeSamFiles.jar INPUT=\$bam_file(s) OUTPUT=\$merged_bam VALIDATION_STRINGENCY=SILENT");

my @header_lines = get_bam_header($split_files[0]);
is_deeply \@header_lines, \@expected_header_lines, 'split bam header is okay';

# check start_from_scratch works correctly
VRPipe::DataElementState->get(pipelinesetup => 2, dataelement => 7)->start_from_scratch();

my @mapping_exists = map { -s $_ ? 1 : 0 } @mapping_files;
is_deeply \@mapping_exists, [1,1,1,1,1], 'mapping files were not deleted after merge element was restarted';

my @split_exists = map { -s $_ ? 1 : 0 } @split_files;
is_deeply \@split_exists, \@split_files_removed, 'correct merge files were removed on start from scratch';

my @release_exists = map { -s $_ ? 1 : 0 } @release_files;
is_deeply \@release_exists, \@release_files_removed, 'correct release files were removed on start from scratch';

ok handle_pipeline(@mapping_files, @split_files, @release_files), 'output files were recreated after a start from scratch';

# if a mapped bam is moved, is everything okay?
my $orig_bam = VRPipe::File->get(path => $mapping_files[4]);
$mapping_pipeline->make_path(dir($orig_bam->dir.'_moved'));
my $moved_bam = VRPipe::File->get(path => file($orig_bam->dir.'_moved', $orig_bam->basename));
$orig_bam->move($moved_bam);
is_deeply [-e $orig_bam->path, -e $moved_bam->path], [undef, 1], 'moved an improved bam';
$mapping_files[4] = file($moved_bam->path);

ok handle_pipeline(@mapping_files, @split_files, @release_files), 'all files are still okay after a bam was moved';

VRPipe::DataElementState->get(pipelinesetup => 2, dataelement => 7)->start_from_scratch();

@mapping_exists = map { -s $_ ? 1 : 0 } @mapping_files;
is_deeply \@mapping_exists, [1,1,1,1,1], 'mapping files were not deleted after merge element was restarted and after bam was moved';

@split_exists = map { -s $_ ? 1 : 0 } @split_files;
is_deeply \@split_exists, \@split_files_removed, 'correct merge files were removed on start from scratch after bam was moved';

@release_exists = map { -s $_ ? 1 : 0 } @release_files;
is_deeply \@release_exists, \@release_files_removed, 'correct release files were removed on start from scratch after bam was moved';

ok handle_pipeline(@mapping_files, @split_files, @release_files), 'output files were recreated from moved bam';

finish;