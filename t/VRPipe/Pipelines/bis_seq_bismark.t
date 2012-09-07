#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES BISMARK_GENOME_FOLDER)],
        required_exe => [qw(fastqc)]
    );
    use TestPipelines;
}
my $output_dir = get_output_dir('bis_seq_bismark-test');
ok my $pipeline = VRPipe::Pipeline->create(name => 'bis_seq_bismark'), 'able to get the bis_seq_bismark pipeline';
my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

is_deeply \@s_names, [qw(fastqc_quality_report trimmomatic bismark sam_sort sam_mark_duplicates bismark_methylation_extractor)], 'the pipeline has the correct steps';

my $pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'bis_seq_bismark_test',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn',
        method  => 'all',
        source  => file(qw(t data fastqc_report_datasource.fofn)),
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        trimmomatic_jar_path  => $ENV{TRIMMOMATIC_JAR_PATH},
        bismark_genome_folder => $ENV{BISMARK_GENOME_FOLDER}
    }
);

my @output_subdirs = output_subdirs(1);

# my $outputfile_1   = file(@output_subdirs, '3_bismark', "2822_6_1.trim", "2822_6_1.trim.fastq_Bismark_mapping_report.txt");
# my $outputfile_2   = file(@output_subdirs, '3_bismark', "2822_6_1.trim", "2822_6_1.trim.fastq_bismark.sam");

my $outputfile_1 = file(@output_subdirs, '6_bismark_methylation_extractor', "CHG_context_2822_6_1.trim.fastq_bismark.sam.sort.markdup.sam.txt");
my $outputfile_2 = file(@output_subdirs, '6_bismark_methylation_extractor', "CpG_context_2822_6_1.trim.fastq_bismark.sam.sort.markdup.sam.txt");
my $outputfile_3 = file(@output_subdirs, '6_bismark_methylation_extractor', "CHH_context_2822_6_1.trim.fastq_bismark.sam.sort.markdup.sam.txt");

my @outputfiles;
push(@outputfiles, $outputfile_1, $outputfile_2, $outputfile_3);
ok handle_pipeline(@outputfiles), 'bismark pipeline ran ok, generating the expected output file';
