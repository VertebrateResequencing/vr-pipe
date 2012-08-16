#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
# NOTE: GSNAP_DB_FOLDER something like mm9 does not require full path
#       gsnap already knows this from installation.

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(gsnap)]
    );
    use TestPipelines;
    use_ok('VRPipe::Steps::gsnap');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('gsnap', 'fastq_files');
is_deeply [$step->id, $step->description], [1, 'Step for GSNAP mapper'], 'GSNAP step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'gsnap',
    datasource => VRPipe::DataSource->create(
        type    => 'delimited',
        method  => 'all_columns',
        options => { delimiter => "\t" },
        source  => file(qw(t data gsnap_datasource_se.fofn))->absolute
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => { gsnap_db => 'mm9', paired_end => 0 }
);

my @output_subdirs = output_subdirs(1);
my $outputfile_1 = file(@output_subdirs, '1_gsnap', 'SRR514151_160_lines.unpaired_uniq');
my @outputfiles;
push(@outputfiles, $outputfile_1);
warn $outputfile_1;
ok handle_pipeline(@outputfiles), 'gsnap pipeline ran ok, generating the expected output file';

#my $testfilecontents   = file( qw(t data ERR032995.concordant_uniq) )->slurp;
#my $outputfilecontents = $outputfile_1->slurp;
#warn $outputfilecontents;
#is($testfilecontents, $outputfilecontents, "output file contain expected data");
finish;
