#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use Path::Class qw(file dir);

BEGIN {
    use Test::Most tests => 10;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPersistentReal;
}

my $scheduler = VRPipe::Scheduler->get();
my $mapping_output_dir = dir($scheduler->output_root, 'pipelines_test_output', 'mapping');
#$scheduler->remove_tree($mapping_output_dir);
#$scheduler->make_path($mapping_output_dir);
my $manager = VRPipe::Manager->get();


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
$scheduler->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);
my $mapping_pipelinesetup = VRPipe::PipelineSetup->get(name => 's_suis mapping',
                                                       datasource => VRPipe::DataSource->get(type => 'sequence_index',
                                                                                             method => 'lane_fastqs',
                                                                                             source => file(qw(t data datasource.sequence_index))),
                                                       output_root => $mapping_output_dir,
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

my @split_fqs = (VRPipe::File->get(path => file($mapping_output_dir, '2822_6', 'fastq_split_se_8000', '2822_6.1.fastq')),
                 VRPipe::File->get(path => file($mapping_output_dir, '2822_6', 'fastq_split_pe_8000', '2822_6_1.1.fastq.gz')),
                 VRPipe::File->get(path => file($mapping_output_dir, '2822_6', 'fastq_split_pe_8000', '2822_6_2.3.fastq.gz')),
                 VRPipe::File->get(path => file($mapping_output_dir, '2823_4', 'fastq_split_pe_8000', '2823_4_2.4.fastq.gz')));

is_deeply [$split_fqs[0]->metadata, $split_fqs[1]->metadata, $split_fqs[2]->metadata, $split_fqs[3]->metadata],
          [{ lane => '2822_6',
             chunk => 0,
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
             mate => file($mapping_output_dir, '2822_6', 'fastq_split_pe_8000', '2822_6_2.1.fastq.gz')->stringify,
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
             mate => file($mapping_output_dir, '2822_6', 'fastq_split_pe_8000', '2822_6_1.3.fastq.gz')->stringify,
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
             mate => file($mapping_output_dir, '2823_4', 'fastq_split_pe_8000', '2823_4_1.4.fastq.gz')->stringify,
             source_fastq => file(qw(t data 2823_4_2.fastq))->absolute->stringify }], 'split files that came out of the fastq_split step have the correct metadata';

my $existing_outputs = 0;
my $existing_sai_outs = 0;
my $existing_sam_outs = 0;
foreach my $lane (qw(2822_6 2822_7 2823_4 8324_8)) {
    if ($lane eq '2822_6') {
        for my $i (1..1) {
            my $fq = file($mapping_output_dir, $lane, 'fastq_split_se_8000', "${lane}.$i.fastq");
            my $sai = $fq.'.sai';
            $existing_outputs += -s $fq ? 1 : 0;
            $existing_sai_outs += -s $sai ? 1 : 0;
            my $sam = file($mapping_output_dir, $lane, 'fastq_split_se_8000', "$lane.sam");
            $existing_sam_outs += -s $sam ? 1 : 0;
        }
        for my $i (1..3) {
            for my $j (1..2) {
                my $fq = file($mapping_output_dir, $lane, 'fastq_split_pe_8000', "${lane}_$j.$i.fastq.gz");
                my $sai = $fq.'.sai';
                $existing_outputs += -s $fq ? 1 : 0;
                $existing_sai_outs += -s $sai ? 1 : 0;
            }
            my $sam = file($mapping_output_dir, $lane, 'fastq_split_pe_8000', "$lane.$i.sam");
            $existing_sam_outs += -s $sam ? 1 : 0;
        }
    }
    else {
        for my $i (1..4) {
            for my $j (1..2) {
                my $fq = file($mapping_output_dir, $lane, 'fastq_split_pe_8000', "${lane}_$j.$i.fastq.gz");
                my $sai = $fq.'.sai';
                $existing_outputs += -s $fq ? 1 : 0;
                $existing_sai_outs += -s $sai ? 1 : 0;
            }
            my $sam = file($mapping_output_dir, $lane, 'fastq_split_pe_8000', "$lane.$i.sam");
            $existing_sam_outs += -s $sam ? 1 : 0;
        }
    }
}

my $recorded_outputs = 0;
foreach my $element_id (1..4) {
    my $fastq_split_outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 2, dataelement => $element_id)->output_files->{split_fastq_files};
    $recorded_outputs += @$fastq_split_outs;
}
is_deeply [$existing_outputs, $recorded_outputs], [31, 31], 'all the split files that should have been created by the fastq_split step exist and were recorded as step outputs';

my $dict_file = VRPipe::File->get(path => file($ref_dir, 'S_suis_P17.fa.dict'));
$existing_outputs = -s $dict_file->path ? 1 : 0;
my %recorded_outputs;
foreach my $element_id (1..4) {
    my $outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 3, dataelement => $element_id)->output_files->{reference_dict};
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
        my $outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 4, dataelement => $element_id)->output_files->{$out_key};
        foreach my $out (@$outs) {
            $recorded_outputs{$out->path->stringify} = 1;
        }
    }
}
$recorded_outputs = keys %recorded_outputs;
is_deeply [$existing_outputs, $recorded_outputs], [8, 8], 'all the ref index files that should have been created by the bwa_index step exist and were recorded as step outputs';

$recorded_outputs = 0;
foreach my $element_id (1..4) {
    my $outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 5, dataelement => $element_id)->output_files->{bwa_sai_files};
    $recorded_outputs += @$outs;
}
is_deeply [$existing_sai_outs, $recorded_outputs], [31, 31], 'all the sai files that should have been created by the bwa_aln_fastq step exist and were recorded as step outputs';

$recorded_outputs = 0;
foreach my $element_id (1..4) {
    my $outs = VRPipe::StepState->get(pipelinesetup => 1, stepmember => 6, dataelement => $element_id)->output_files->{bwa_sam_files};
    $recorded_outputs += @$outs;
}
is_deeply [$existing_sam_outs, $recorded_outputs], [16, 16], 'all the sam files that should have been created by the bwa_sam step exist and were recorded as step outputs';

done_testing;
exit;

sub handle_pipeline {
    my $give_up = 200;
    while (! $manager->trigger) {
        last if $give_up-- <= 0;
        $manager->handle_submissions;
        sleep(1);
    }
    my $all_created = 1;
    foreach my $ofile (@_) {
        unless (-s $ofile) {
            warn "$ofile is missing\n";
            $all_created = 0;
        }
    }
    return $all_created;
}