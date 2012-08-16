#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)], # require bismark path
        required_exe => [qw(methylation_extractor)]
    );
    use TestPipelines;
    use_ok('VRPipe::Steps::bismark_methylation_extractor');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bismark_methylation_extractor', 'sam_file');
is_deeply [$step->id, $step->description], [1, 'Step for the Methylation Extractor tool bundled with Bismark to produce a table of methylation calls by position'], 'Bismark methylation extractor step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'bismark_methylation_extractor',
    datasource => VRPipe::DataSource->create(
        type    => 'delimited',
        method  => 'all_columns',
        options => { delimiter => "\t" },
        source  => file(qw(t data bismark_meth_exr_datasource.fofn))->absolute
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => { paired_end => 1 }
);

my @output_subdirs = output_subdirs(1);
my $outputfile_1   = file(@output_subdirs, '1_bismark_methylation_extractor', "CHG_context_2822_6_1.fastq_bismark_pe.sam.txt");
my $outputfile_2   = file(@output_subdirs, '1_bismark_methylation_extractor', "CpG_context_2822_6_1.fastq_bismark_pe.sam.txt");
my $outputfile_3   = file(@output_subdirs, '1_bismark_methylation_extractor', "CHH_context_2822_6_1.fastq_bismark_pe.sam.txt");
my @outputfiles;
push @outputfiles, $outputfile_1, $outputfile_2, $outputfile_3;
ok handle_pipeline(@outputfiles), 'bismark methylation extractor pipeline ran ok, generating the expected file';
