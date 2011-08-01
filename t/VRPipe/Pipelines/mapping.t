#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 14;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPipelines;
}

my $mapping_output_dir = get_output_dir('mapping');

ok my $mapping_pipeline = VRPipe::Pipeline->get(name => 'fastq_mapping_with_bwa'), 'able to get a pre-written pipeline';
#TODO: {
#    local $TODO = "not all steps implemented yet";
#    my @s_names;
#    foreach my $stepmember ($mapping_pipeline->steps) {
#        push(@s_names, $stepmember->step->name);
#    }
#    is_deeply \@s_names, [qw(fastq_metadata fastq_split bwa_index bwa_aln bwa_sam sam_to_fixed_bam bam_merge_lane_splits bam_stats)], 'the pipeline has the correct steps';
#    
#    $mapping_pipeline = VRPipe::Pipeline->get(name => 'mapping');
#    @s_names = ();
#    foreach my $stepmember ($mapping_pipeline->steps) {
#        push(@s_names, $stepmember->step->name);
#    }
#    is_deeply \@s_names, [qw(fastq_metadata fastq_split bwa_index bwa_aln bwa_sam sam_to_fixed_bam bam_merge_lane_splits bam_stats)], 'the pipeline has the correct steps after a second retrieval';
#}

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
                                                                   bwa_index_cmd => 'bwa index -a is',
                                                                   cleanup => 0});


# make a second setup in a different dir; whilst we don't have any tests
# specific to this setup, creating it is sufficient to reveal problems with
# forking in Manager->trigger
my $mapping_output_dir2 = get_output_dir('mapping_cleanup');
VRPipe::PipelineSetup->get(name => 's_suis mapping',
                           datasource => VRPipe::DataSource->get(type => 'sequence_index',
                                                                 method => 'lane_fastqs',
                                                                 source => file(qw(t data datasource.sequence_index)),
                                                                 options => { local_root_dir => dir(".")->absolute->stringify }),
                           output_root => $mapping_output_dir2,
                           pipeline => $mapping_pipeline,
                           options => {fastq_chunk_size => 8000,
                                       reference_fasta => $ref_fa,
                                       reference_assembly_name => 'SSuis1',
                                       reference_public_url => 'ftp://s.suis.com/ref.fa',
                                       reference_species => 'S.Suis',
                                       bwa_index_cmd => 'bwa index -a is'});

my @mapping_output_files;

handle_pipeline(@mapping_output_files);

is_deeply [VRPipe::File->get(path => file(qw(t data 8324_8_1.fastq))->absolute)->md5,
           VRPipe::File->get(path => file(qw(t data 8324_8_1.fastq))->absolute)->metadata,
           VRPipe::File->get(path => file(qw(t data 8324_8_2.fastq))->absolute)->metadata],
          ['b90eb321fb12b47800132b237a6a3722',
           { lane => '8324_8',
             study => 'STUDY01',
             study_name => 'my study name',
             center_name => 'SC',
             sample_id => 'SAMPLEID02',
             sample => 'SAMPLE02',
             population => 'POP',
             platform => 'ILLUMINA',
             library => 'LIB03',
             insert_size => 200,
             withdrawn => 0,
             reads => 250,
             bases => 15250,
             avg_read_length => '61.00',
             analysis_group => 'low coverage',
             paired => 1,
             mate => file(qw(t data 8324_8_2.fastq))->absolute->stringify },
           { expected_md5 => 'c8ecf1168b710adde2a3b45f98849780',
             lane => '8324_8',
             study => 'STUDY01',
             study_name => 'my study name',
             center_name => 'SC',
             sample_id => 'SAMPLEID02',
             sample => 'SAMPLE02',
             population => 'POP',
             platform => 'ILLUMINA',
             library => 'LIB03',
             insert_size => 200,
             withdrawn => 0,
             reads => 250,
             bases => 13500,
             avg_read_length => '54.00',
             analysis_group => 'low coverage',
             paired => 2,
             mate => file(qw(t data 8324_8_1.fastq))->absolute->stringify }], 'fastqs that went through the first step of the mapping pipeline have the correct metadata';

my @split_fqs = (VRPipe::File->get(path => file($mapping_output_dir, 1, 'fastq_split', 'se_8000', '2822_6.1.fastq')),
                 VRPipe::File->get(path => file($mapping_output_dir, 1, 'fastq_split', 'pe_8000', '2822_6_1.1.fastq.gz')),
                 VRPipe::File->get(path => file($mapping_output_dir, 1, 'fastq_split', 'pe_8000', '2822_6_2.3.fastq.gz')),
                 VRPipe::File->get(path => file($mapping_output_dir, 3, 'fastq_split', 'pe_8000', '2823_4_2.4.fastq.gz')));

is_deeply [$split_fqs[0]->metadata, $split_fqs[1]->metadata, $split_fqs[2]->metadata, $split_fqs[3]->metadata],
          [{ lane => '2822_6',
             chunk => 1,
             study => 'STUDY01',
             study_name => 'my study name',
             center_name => 'SC',
             sample_id => 'SAMPLEID01',
             sample => 'SAMPLE01',
             population => 'POP',
             platform => 'ILLUMINA',
             library => 'LIB01',
             insert_size => 200,
             withdrawn => 0,
             reads => 50,
             bases => 3050,
             avg_read_length => '61.00',
             analysis_group => 'low coverage',
             paired => 0,
             source_fastq => file(qw(t data 2822_6.fastq))->absolute->stringify },
           { lane => '2822_6',
             chunk => 1,
             study => 'STUDY01',
             study_name => 'my study name',
             center_name => 'SC',
             sample_id => 'SAMPLEID01',
             sample => 'SAMPLE01',
             population => 'POP',
             platform => 'ILLUMINA',
             library => 'LIB01',
             insert_size => 200,
             withdrawn => 0,
             reads => 69,
             bases => 4209,
             avg_read_length => '61.00',
             analysis_group => 'low coverage',
             paired => 1,
             mate => file($mapping_output_dir, 1, 'fastq_split', 'pe_8000', '2822_6_2.1.fastq.gz')->stringify,
             source_fastq => file(qw(t data 2822_6_1.fastq))->absolute->stringify },
           { lane => '2822_6',
             chunk => 3,
             study => 'STUDY01',
             study_name => 'my study name',
             center_name => 'SC',
             sample_id => 'SAMPLEID01',
             sample => 'SAMPLE01',
             population => 'POP',
             platform => 'ILLUMINA',
             library => 'LIB01',
             insert_size => 200,
             withdrawn => 0,
             reads => 62,
             bases => 3348,
             avg_read_length => '54.00',
             analysis_group => 'low coverage',
             paired => 2,
             mate => file($mapping_output_dir, 1, 'fastq_split', 'pe_8000', '2822_6_1.3.fastq.gz')->stringify,
             source_fastq => file(qw(t data 2822_6_2.fastq))->absolute->stringify },
           { lane => '2823_4',
             chunk => 4,
             study => 'STUDY01',
             study_name => 'my study name',
             center_name => 'SC',
             sample_id => 'SAMPLEID01',
             sample => 'SAMPLE01',
             population => 'POP',
             platform => 'ILLUMINA',
             library => 'LIB02',
             insert_size => 200,
             withdrawn => 0,
             reads => 43,
             bases => 2322,
             avg_read_length => '54.00',
             analysis_group => 'low coverage',
             paired => 2,
             mate => file($mapping_output_dir, 3, 'fastq_split', 'pe_8000', '2823_4_1.4.fastq.gz')->stringify,
             source_fastq => file(qw(t data 2823_4_2.fastq))->absolute->stringify }], 'split files that came out of the fastq_split step have the correct metadata';

my $existing_outputs = 0;
my $existing_sai_outs = 0;
my $existing_sam_outs = 0;
my $existing_split_bam_outs = 0;
my $existing_final_bam_outs = 0;
my @final_bams;
my $element_num = 0;
foreach my $lane (qw(2822_6 2822_7 2823_4 8324_8)) {
    $element_num++;
    if ($lane eq '2822_6') {
        for my $i (1..1) {
            my $fq = file($mapping_output_dir, $element_num, 'fastq_split', 'se_8000', "${lane}.$i.fastq");
            my $sai = $fq.'.sai';
            $existing_outputs += -s $fq ? 1 : 0;
            $existing_sai_outs += -s $sai ? 1 : 0;
            my $sam = file($mapping_output_dir, $element_num, 'bwa_sam', "$lane.se.1.sam");
            $existing_sam_outs += -s $sam ? 1 : 0;
            my $bam = file($mapping_output_dir, $element_num, 'sam_to_fixed_bam', "$lane.pe.$i.bam");
            $existing_split_bam_outs += -s $bam ? 1 : 0;
        }
        for my $i (1..3) {
            for my $j (1..2) {
                my $fq = file($mapping_output_dir, $element_num, 'fastq_split', 'pe_8000', "${lane}_$j.$i.fastq.gz");
                my $sai = $fq.'.sai';
                $existing_outputs += -s $fq ? 1 : 0;
                $existing_sai_outs += -s $sai ? 1 : 0;
            }
            my $sam = file($mapping_output_dir, $element_num, 'bwa_sam', "$lane.pe.$i.sam");
            $existing_sam_outs += -s $sam ? 1 : 0;
            my $bam = file($mapping_output_dir, $element_num, 'sam_to_fixed_bam', "$lane.pe.$i.bam");
            $existing_split_bam_outs += -s $bam ? 1 : 0;
        }
        
        foreach my $ended ('se', 'pe') {
            my $bam = file($mapping_output_dir, $element_num, 'bam_merge_lane_splits', "$lane.$ended.bam");
            $existing_final_bam_outs += -s $bam ? 1 : 0;
            push(@final_bams, VRPipe::File->get(path => $bam));
        }
    }
    else {
        for my $i (1..4) {
            for my $j (1..2) {
                my $fq = file($mapping_output_dir, $element_num, 'fastq_split', 'pe_8000', "${lane}_$j.$i.fastq.gz");
                my $sai = $fq.'.sai';
                $existing_outputs += -s $fq ? 1 : 0;
                $existing_sai_outs += -s $sai ? 1 : 0;
            }
            my $sam = file($mapping_output_dir, $element_num, 'bwa_sam', "$lane.pe.$i.sam");
            $existing_sam_outs += -s $sam ? 1 : 0;
            my $bam = file($mapping_output_dir, $element_num, 'sam_to_fixed_bam', "$lane.pe.$i.bam");
            $existing_split_bam_outs += -s $bam ? 1 : 0;
        }
        
        my $bam = file($mapping_output_dir, $element_num, 'bam_merge_lane_splits', "$lane.pe.bam");
        $existing_final_bam_outs += -s $bam ? 1 : 0;
        push(@final_bams, VRPipe::File->get(path => $bam));
    }
}

my $recorded_outputs = 0;
foreach my $element_id (1..4) {
    my $fastq_split_outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 4, dataelement => $element_id)->output_files->{split_fastq_files};
    $recorded_outputs += @$fastq_split_outs;
}
is_deeply [$existing_outputs, $recorded_outputs], [31, 31], 'all the split files that should have been created by the fastq_split step exist and were recorded as step outputs';

my $dict_file = VRPipe::File->get(path => file($ref_dir, 'S_suis_P17.fa.dict'));
$existing_outputs = -s $dict_file->path ? 1 : 0;
my %recorded_outputs;
foreach my $element_id (1..4) {
    my $outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 1, dataelement => $element_id)->output_files->{reference_dict};
    foreach my $out (@$outs) {
        $recorded_outputs{$out->path->stringify}++;
    }
}
$recorded_outputs = keys %recorded_outputs;
is_deeply [$existing_outputs, $recorded_outputs, $recorded_outputs{$dict_file->path}], [1, 1, 4], 'the .dict file from sequence_dictionary step was made and recorded as the step output';
is $dict_file->slurp, "\@HD\tVN:1.0\tSO:unsorted
\@SQ\tSN:fake_chr1\tLN:290640\tM5:55f9584cf1f4194f13bbdc0167e0a05f\tUR:ftp://s.suis.com/ref.fa\tAS:SSuis1\tSP:S.Suis
\@SQ\tSN:fake_chr2\tLN:1716851\tM5:6dd2836053e5c4bd14ad49b5b2f2eb88\tUR:ftp://s.suis.com/ref.fa\tAS:SSuis1\tSP:S.Suis\n", '.dict file content is correct';

$existing_outputs = 0;
foreach my $suffix (qw(amb ann bwt pac rbwt rpac rsa sa)) {
    $existing_outputs += -s file($ref_dir, 'S_suis_P17.fa.'.$suffix) ? 1 : 0;
}
%recorded_outputs = ();
foreach my $element_id (1..4) {
    foreach my $out_key (qw(bwa_index_binary_files bwa_index_text_files)) {
        my $outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 2, dataelement => $element_id)->output_files->{$out_key};
        foreach my $out (@$outs) {
            $recorded_outputs{$out->path->stringify} = 1;
        }
    }
}
$recorded_outputs = keys %recorded_outputs;
is_deeply [$existing_outputs, $recorded_outputs], [8, 8], 'all the ref index files that should have been created by the bwa_index step exist and were recorded as step outputs';

$recorded_outputs = 0;
foreach my $element_id (1..4) {
    my $outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 6, dataelement => $element_id)->output_files->{bwa_sai_files};
    $recorded_outputs += @$outs;
}
is_deeply [$existing_sai_outs, $recorded_outputs], [31, 31], 'all the sai files that should have been created by the bwa_aln_fastq step exist and were recorded as step outputs';

$recorded_outputs = 0;
foreach my $element_id (1..4) {
    my $outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 7, dataelement => $element_id)->output_files->{bwa_sam_files};
    $recorded_outputs += @$outs;
}
is_deeply [$existing_sam_outs, $recorded_outputs], [16, 16], 'all the sam files that should have been created by the bwa_sam step exist and were recorded as step outputs';

$recorded_outputs = 0;
foreach my $element_id (1..4) {
    my $outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 8, dataelement => $element_id)->output_files->{fixed_bam_files};
    $recorded_outputs += @$outs;
}
is_deeply [$existing_split_bam_outs, $recorded_outputs], [16, 16], 'all the bam files that should have been created by the sam_to_fixed_bam step exist and were recorded as step outputs';

$recorded_outputs = 0;
foreach my $element_id (1..4) {
    my $outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 9, dataelement => $element_id)->output_files->{merged_lane_bams};
    $recorded_outputs += @$outs;
}
is_deeply [$existing_final_bam_outs, $recorded_outputs], [5, 5], 'all the bam files that should have been created by the bam_merge_lane_splits step exist and were recorded as step outputs';

my @final_bam_metas = map { $_->metadata } sort { $a->path cmp $b->path } @final_bams;
is_deeply [@final_bam_metas],
          [{ lane => '2822_6',
             study => 'STUDY01',
             center_name => 'SC',
             sample => 'SAMPLE01',
             platform => 'ILLUMINA',
             library => 'LIB01',
             insert_size => 200,
             reads => 400,
             bases => 23000,
             paired => 1,
             mapped_fastqs => join(',', file(qw(t data 2822_6_1.fastq))->absolute->stringify, file(qw(t data 2822_6_2.fastq))->absolute->stringify) },
           { lane => '2822_6',
             study => 'STUDY01',
             center_name => 'SC',
             sample => 'SAMPLE01',
             platform => 'ILLUMINA',
             library => 'LIB01',
             insert_size => 200,
             reads => 50,
             bases => 3050,
             paired => 0,
             mapped_fastqs => file(qw(t data 2822_6.fastq))->absolute->stringify },
           { lane => '2822_7',
             study => 'STUDY01',
             center_name => 'SC',
             sample => 'SAMPLE01',
             platform => 'ILLUMINA',
             library => 'LIB01',
             insert_size => 200,
             reads => 500,
             bases => 28750,
             paired => 1,
             mapped_fastqs => join(',', file(qw(t data 2822_7_1.fastq))->absolute->stringify, file(qw(t data 2822_7_2.fastq))->absolute->stringify) },
           { lane => '2823_4',
             study => 'STUDY01',
             center_name => 'SC',
             sample => 'SAMPLE01',
             platform => 'ILLUMINA',
             library => 'LIB02',
             insert_size => 200,
             reads => 500,
             bases => 28750,
             paired => 1,
             mapped_fastqs => join(',', file(qw(t data 2823_4_1.fastq))->absolute->stringify, file(qw(t data 2823_4_2.fastq))->absolute->stringify) },
           { lane => '8324_8',
             study => 'STUDY01',
             center_name => 'SC',
             sample => 'SAMPLE02',
             platform => 'ILLUMINA',
             library => 'LIB03',
             insert_size => 200,
             reads => 500,
             bases => 28750,
             paired => 1,
             mapped_fastqs => join(',', file(qw(t data 8324_8_1.fastq))->absolute->stringify, file(qw(t data 8324_8_2.fastq))->absolute->stringify) }], 'final bam files have the correct metadata';

is_deeply [VRPipe::StepState->get(pipelinesetup => 1, stepmember => 2, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 6, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 7, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 8, dataelement => 1)->cmd_summary->summary],
          ['bwa index -a is $reference_fasta',
           'bwa aln -q 15 -f $sai_file $reference_fasta $fastq_file',
           'bwa sampe -a 600 -r $rg_line -f $sam_file $reference_fasta $sai_file(s) $fastq_file(s)',
           'samtools view -bSu $sam_file | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - $reference_fasta > $fixed_bam_file'],
          'cmd summaries for the major steps were as expected';

done_testing;
exit;