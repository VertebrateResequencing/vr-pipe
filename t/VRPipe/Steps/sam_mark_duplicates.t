#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 2;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)] # require bismark path
    );
    use TestPipelines;
    use_ok('VRPipe::Steps::sam_mark_duplicates');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('sam_mark_duplicates', 'sam_files');
is_deeply [$step->id, $step->description], [1, 'Mark duplicates in a sam files using picard'], 'Bismark methylation extractor step created and has correct description';



#my $setup = VRPipe::PipelineSetup->create(name       => 'sam_mark_duplicates',
#                                          datasource => VRPipe::DataSource->create(type    => 'delimited',
#                                                                                   method  => 'all_columns',
#                                                                                   options => { delimiter => "\t" },
#                                                                                   source  => file(qw(t data bismark_meth_exr_datasource.fofn))->absolute),
#                                          output_root => $output_dir,
#                                          pipeline    => $pipeline,
#                                          options     => { paired_end => 1 });
#my @output_subdirs = output_subdirs(1);
#my $outputfile_1   = file(@output_subdirs, '1_sam_mark_duplicates', "CHG_context_2822_6_1.fastq_bismark_pe.sam.txt");
#my @outputfiles;
#push @outputfiles, $outputfile_1;
#warn $outputfile_1;
#ok handle_pipeline(@outputfiles), 'bismark methylation extractor pipeline ran ok, generating the expected file';
