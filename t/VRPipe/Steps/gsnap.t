#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
# NOTE: GSNAP_DB_FOLDER something like mm9 does not require full path
#       gsnap already knows this from installation.

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(gsnap)]
    );
    use TestPipelines;
    use_ok('VRPipe::Steps::gsnap');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('gsnap', 'fastq_files');
my $gmap_build_dir = dir(qw(t data sacCer3 gmap_build));
SKIP: {
    skip "no genome directory found to run test, try running gmap_build test", 2 if (!$gmap_build_dir->stat);
    is_deeply [$step->id, $step->description], [1, 'Step for GSNAP mapper'], 'GSNAP step created and has correct description';
    
    my $setup = VRPipe::PipelineSetup->create(
        name       => 'gsnap',
        datasource => VRPipe::DataSource->create(
            type    => 'delimited',
            method  => 'all_columns',
            options => { delimiter => "\t" },
            source  => file(qw(t data gsnap_datasource.fofn))->absolute
        ),
        output_root => dir($output_dir)->absolute->stringify,
        pipeline    => $pipeline,
        options     => { gsnap_db => 'mygenome', paired_end => 1, gsnap_genome_dir => $gmap_build_dir->absolute->stringify }
    );
    
    my @output_subdirs = output_subdirs(1);
    my $outputfile_1 = file(@output_subdirs, '1_gsnap', 'SRR364359_1_100000lines.concordant_uniq');
    my @outputfiles;
    push(@outputfiles, $outputfile_1);
    ok handle_pipeline(@outputfiles), 'gsnap step ran ok, generating the expected output file';
    
    my $testfilecontents = file(qw(t data SRR364359_1_100000lines.concordant_uniq));
    
    #the order of records may not be the same from one run to the next so
    #need to sort them
    my @testfilecontents   = sort $testfilecontents->slurp;
    my @outputfilecontents = sort $outputfile_1->slurp;
    is_deeply(\@testfilecontents, \@outputfilecontents, "output file contain expected data");

}
finish;
