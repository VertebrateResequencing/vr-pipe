#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 6;
    # this test is Sanger-specific, only the author needs to run it
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB)],
                    required_exe => [qw(iget iquest)]);
    use TestPipelines;
    
    use_ok('VertRes::Utils::VRTrackFactory');
}

# setup a little VRTrack db that has its files in irods
my %cd = VertRes::Utils::VRTrackFactory->connection_details('rw');
open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
my @sql = VRPipe::File->get(path => file(qw(t data vrtrack_yeast_irods.sql))->absolute)->slurp;
foreach my $sql (@sql) {
    print $mysqlfh $sql;
}
close($mysqlfh);

my $output_dir = get_output_dir('bam_import_from_irods_and_vrtrack_qc_wgs');
my $irods_dir = dir($output_dir, 'irods_import')->stringify;

ok my $ds = VRPipe::DataSource->get(type => 'vrtrack',
                                    method => 'lane_bams',
                                    source => $ENV{VRPIPE_VRTRACK_TESTDB},
                                    options => {local_root_dir => $irods_dir}), 'could create a vrtrack datasource';
my $results = 0;
foreach my $element (@{$ds->elements}) {
    $results++;
}
is $results, 29, 'got correct number of bams from the vrtrack db';

ok my $import_qc_pipeline = VRPipe::Pipeline->get(name => 'bam_import_from_irods_and_vrtrack_qc_wgs'), 'able to get the bam_import_from_irods_and_vrtrack_qc_wgs pipeline';
my @s_names;
foreach my $stepmember ($import_qc_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(irods_get_files_by_basename bamcheck plot_bamcheck bamcheck_rmdup bamcheck_stats_output vrtrack_update_mapstats)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data pombe_ref.fa));
my $ref_dir = dir($output_dir, 'ref');
$import_qc_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'pombe_ref.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $stats_fa_source = $ref_fa_source.'.stats';
my $ref_fa_stats = $ref_fa.'.stats';
copy($stats_fa_source, $ref_fa_stats);

VRPipe::PipelineSetup->get(name => 'pombe import and qc',
                           datasource => $ds,
                           output_root => $output_dir,
                           pipeline => $import_qc_pipeline,
                           options => {reference_fasta => $ref_fa,
                                       reference_assembly_name => 'SPombe1',
                                       reference_public_url => 'ftp://s.pombe.com/ref.fa',
                                       reference_species => 'S.Pombe',
                                       reference_fasta_stats => $ref_fa_stats,
                                       bamcheck_options => '-q 20 -r', # -q20 -d for bamcheck_rmdup_options?
                                       vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB},
                                       cleanup => 0});

my @irods_files;
foreach my $num (qw(1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26 27 28 29 30 31)) {
    push(@irods_files, file($irods_dir, '7369_5#'.$num.'.bam'));
}
ok handle_pipeline(@irods_files), 'pipeline ran ok';


# we'll also test the whole chain of pipelines we typically run in
# VertebrateResequencing at the Sanger: improvement followed by genotype check
# followed by merge up to the sample level

#...

finish;