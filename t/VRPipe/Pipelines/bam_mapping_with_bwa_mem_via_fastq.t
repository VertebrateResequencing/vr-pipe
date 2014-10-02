#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 8;
    use VRPipeTest (
        required_env => 'VRPIPE_TEST_PIPELINES',
        required_exe => [qw(bwa samtools bamtofastq split)]
    );
    use TestPipelines;
}

## first bwa mem pipeline
my $mapping_output_dir = get_output_dir('bam_mapping_with_bwa_mem_via_fastq');

ok my $mapping_pipeline = VRPipe::Pipeline->create(name => 'bam_mapping_with_bwa_mem_via_fastq'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($mapping_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(fasta_index sequence_dictionary bwa_index bam_metadata bam_shuffle_by_name bam_to_fastq fastq_split bwa_mem_fastq sam_to_fixed_bam bam_merge_lane_splits bam_index)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($mapping_output_dir, 'ref');
$mapping_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);
my $mapping_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'bam_mapping_with_bwa_mem_via_fastq',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data datasource.bam_fofn))->absolute
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
        bwa_mem_fastq_memory          => 150,
        bwa_mem_fastq_time            => 1,
        sam_to_fixed_bam_memory       => 150,
        sam_to_fixed_bam_time         => 1,
        bam_merge_lane_splits_memory  => 150,
        bam_merge_lane_splits_time    => 1,
        cleanup                       => 0
    }
);
ok handle_pipeline(), 'pipeline ran ok';

## second bwa mem pipeline
my $mapping_output_dir2 = get_output_dir('bam_remapping_with_bwa_mem');

ok my $mapping_pipeline2 = VRPipe::Pipeline->create(name => 'bam_remapping_with_bwa_mem'), 'able to get a pre-written pipeline';

my @s_names2;
foreach my $stepmember ($mapping_pipeline2->step_members) {
    push(@s_names2, $stepmember->step->name);
}
is_deeply \@s_names2, [qw(fasta_index sequence_dictionary bwa_index bam_metadata bamtofastq bwa_mem_to_bam bam_merge_lane_splits bam_index)], 'the pipeline has the correct steps';

my $mapping_pipelinesetup2 = VRPipe::PipelineSetup->create(
    name       => 'bam_remapping_with_bwa_mem',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data datasource.bam_fofn))->absolute
    ),
    output_root => $mapping_output_dir2,
    pipeline    => $mapping_pipeline2,
    options     => {
        reference_fasta              => $ref_fa,
        reference_assembly_name      => 'SSuis1',
        reference_public_url         => 'ftp://s.suis.com/ref.fa',
        reference_species            => 'S.Suis',
        bwa_index_options            => '-a is',
        fastq_chunk_size             => 8000,
        sequence_dictionary_memory   => 150,
        sequence_dictionary_time     => 1,
        bwa_mem_memory               => 150,
        bwa_mem_time                 => 1,
        bam_merge_lane_splits_memory => 150,
        bam_merge_lane_splits_time   => 1,
        bwa_index_memory             => 150,
        bwa_index_time               => 1,
        cleanup                      => 0
    }
);
ok handle_pipeline(), 'pipeline ran ok';

### checking metadata of final bam files for each setup. The pipelines should
## produce the expected metadata if they run correctly.
for my $setup_num (1, 2) {
    my $bam_step_num = 10; # bam lanes merge step #
    $bam_step_num = 7 unless ($setup_num == 1);
    my $bamtofastq_dir = "6_bam_to_fastq"; #used for setup 1 only
    my $shuf           = ".shuf";
    $shuf = '' unless ($setup_num == 1);
    
    my %mapped_fastqs;
    my %source_bam;
    my @final_bams;
    my $source_type;
    my %source_values;
    my $element_num = 0;
    foreach my $lane (qw(2822_6 2822_7 2823_4 8324_8)) {
        if ($lane eq '2822_6') {
            foreach my $ended ('se', 'pe') {
                $element_num++;
                my @output_subdirs = output_subdirs($element_num, $setup_num);
                my $bam = file(@output_subdirs, "${bam_step_num}_bam_merge_lane_splits", "$lane.$ended.bam");
                push(@final_bams, VRPipe::File->create(path => $bam));
                my @fqs = ();
                if ($ended eq 'se') {
                    push(@fqs, file(@output_subdirs, $bamtofastq_dir, "$lane.$ended${shuf}_M.fastq"));
                }
                else {
                    for my $j (1 .. 2) {
                        push(@fqs, file(@output_subdirs, $bamtofastq_dir, "$lane.$ended${shuf}_$j.fastq"));
                    }
                }
                $mapped_fastqs{"$lane.$ended"} = join(',', @fqs);
                $source_bam{"$lane.$ended"} = file("t", "data", "remapping_bams", "$lane.$ended.bam")->absolute;
            }
        }
        else {
            $element_num++;
            my @output_subdirs = output_subdirs($element_num, $setup_num);
            my $bam = file(@output_subdirs, "${bam_step_num}_bam_merge_lane_splits", "$lane.pe.bam");
            push(@final_bams, VRPipe::File->create(path => $bam));
            my @fqs = ();
            for my $j (1 .. 2) {
                push(@fqs, file(@output_subdirs, $bamtofastq_dir, "$lane.pe${shuf}_$j.fastq"));
            }
            $mapped_fastqs{$lane} = join(',', @fqs);
            $source_bam{$lane} = file("t", "data", "remapping_bams", "$lane.pe.bam")->absolute;
        }
    }
    if ($setup_num == 1) {
        $source_type   = "mapped_fastqs";
        %source_values = %mapped_fastqs;
    }
    else {
        $source_type   = "source_bam";
        %source_values = %source_bam;
    }
    
    my @final_bam_metas = map { $_->metadata } @final_bams;
    is_deeply [@final_bam_metas],
      [{
            bases        => 3050,
            center_name  => 'SC',
            lane         => '2822_6',
            library      => 'LIB01',
            $source_type => $source_values{'2822_6.se'},
            paired       => 0,
            platform     => 'ILLUMINA',
            reads        => 50,
            sample       => 'SAMPLE01',
            study        => 'STUDY01'
        },
        {
            bases        => 23000,
            center_name  => 'SC',
            lane         => '2822_6',
            library      => 'LIB01',
            $source_type => $source_values{'2822_6.pe'},
            paired       => 1,
            platform     => 'ILLUMINA',
            reads        => 400,
            sample       => 'SAMPLE01',
            study        => 'STUDY01'
        },
        {
            bases        => 28750,
            center_name  => 'SC',
            lane         => '2822_7',
            library      => 'LIB01',
            $source_type => $source_values{'2822_7'},
            paired       => 1,
            platform     => 'ILLUMINA',
            reads        => 500,
            sample       => 'SAMPLE01',
            study        => 'STUDY01'
        },
        {
            bases        => 28750,
            center_name  => 'SC',
            lane         => '2823_4',
            library      => 'LIB02',
            $source_type => $source_values{'2823_4'},
            paired       => 1,
            platform     => 'ILLUMINA',
            reads        => 500,
            sample       => 'SAMPLE01',
            study        => 'STUDY01'
        },
        {
            bases        => 28750,
            center_name  => 'SC',
            lane         => '8324_8',
            library      => 'LIB03',
            $source_type => $source_values{'8324_8'},
            paired       => 1,
            platform     => 'ILLUMINA',
            reads        => 500,
            sample       => 'SAMPLE02',
            study        => 'STUDY01'
        }
      ],
      "final bam files for setup $setup_num have the correct metadata";
}

finish;
