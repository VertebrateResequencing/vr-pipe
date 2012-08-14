#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)] # require bismark path
    );
    use TestPipelines;
    use_ok('VRPipe::Steps::sam_sort');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('sam_sort', 'sam_file');
is_deeply [$step->id, $step->description], [1, "Sort a sam file using - grep '^\@' file_to_sort.sam > sorted_file.sam; grep-v '^\@' file_to_sort.sam | sort -k 3,3 -k4,4n >> sorted_file.sam"], 'Sort sam step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(name       => 'sam_sort',
                                          datasource => VRPipe::DataSource->create(type    => 'delimited',
                                                                                   method  => 'all_columns',
                                                                                   options => { delimiter => "\t" },
                                                                                   source  => file(qw(t data bismark_meth_exr_datasource.fofn))->absolute),
                                          output_root => $output_dir,
                                          pipeline    => $pipeline,
                                          options     => {});
my @output_subdirs = output_subdirs(1);
my $outputfile_1 = file(@output_subdirs, '1_sam_sort', "2822_6_1.fastq_bismark_pe.sam.sort.sam");

my @outputfiles;
push @outputfiles, $outputfile_1;
ok handle_pipeline(@outputfiles), 'bismark methylation extractor pipeline ran ok, generating the expected file';

# is it sorted
my $testfilecontents   = file(qw( t data 2822_6_1.fastq_bismark_pe.sorted.sam ))->slurp;
my $outputfilecontents = $outputfile_1->slurp;
is($testfilecontents, $outputfilecontents, 'file is correctly sorted by co-ordinate');
