#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 7;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB)],
        required_exe => [qw(iget iquest)]
    );
    use TestPipelines;
    
    use_ok('VRTrack::Factory');
}

my %cd = VRTrack::Factory->connection_details('rw');
open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
my @sql = VRPipe::File->create(path => file(qw(t data vrtrack_hipsci_qc1_pilot.sql))->absolute)->slurp;
foreach my $sql (@sql) {
    print $mysqlfh $sql;
}
close($mysqlfh);

# setup pipeline
my $output_dir = get_output_dir('genome_studio_import_genotype_files');
my $irods_dir = dir($output_dir, 'irods_import')->stringify;

#setup vrtrack datasource
ok my $ds = VRPipe::DataSource->create(
    type    => 'vrtrack',
    method  => 'analysis_gtc',
    source  => $ENV{VRPIPE_VRTRACK_TESTDB},
    options => { local_root_dir => $irods_dir }
  ),
  'could create a vrtrack datasource';

#check correct number of gtc file retrieved
my $results = 0;
foreach my $element (@{ get_elements($ds) }) {
    $results++;
}
is $results, 5, 'got correct number of gtc files from the vrtrack db';

#check pieline has correct steps
ok my $pipeline = VRPipe::Pipeline->create(name => 'genome_studio_import_genotype_files'), 'able to get the genome_studio_import_genotype_files pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(irods_get_files_by_basename split_genome_studio_genotype_files)], 'the pipeline has the correct steps';

#create external genotype gzip file for testing to override the path in gtc file metadata
my $external_gzip_source = file(qw(t data test_pilot_genotyping.fcr.txt.gz));
my $gzip_dir = dir($output_dir, 'external_gzip');
$pipeline->make_path($gzip_dir);
my $external_gzip_file = file($gzip_dir, 'test_pilot_genotyping.fcr.txt.gz')->stringify;
copy($external_gzip_source, $external_gzip_file);

# create pipeline setup
VRPipe::PipelineSetup->create(
    name        => 'gtc import and qc',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        vrtrack_db         => $ENV{VRPIPE_VRTRACK_TESTDB},
        irods_get_zone     => 'archive',
        external_gzip_file => $external_gzip_file,
        cleanup            => 1
    }
);

#get arrays of output files
my @irods_files;
my @lanes = qw(9274735028_R02C01 9274735028_R03C01 9274735028_R04C01 9274735028_R05C01 9274735028_R06C01);
foreach my $lane (@lanes) {
    push(@irods_files, file($irods_dir, $lane . '.gtc'));
}

my @genotype_files;
my $element_id = 0;
foreach my $sample (qw(271298_B03_hipscigt5466711 271298_A03_hipscigt5466710 271298_H02_hipscigt5466709 271298_G02_hipscigt5466708 271298_F02_hipscigt5466707)) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@genotype_files, file(@output_subdirs, '2_split_genome_studio_genotype_files', $sample . '.genotyping.fcr.txt'));
}

#run pipeline and check outputs
ok handle_pipeline(@genotype_files), 'bam import from irods and split genome studio genotype file pipeline ran ok';

#check genotype file metadata
my $meta = VRPipe::File->get(path => $genotype_files[0])->metadata;
is_deeply $meta,
  {
    'analysis_uuid' => 'b37c6db3-c2bb-40db-aa1b-48340681ebe1',
    'bases'         => '0',
    'withdrawn'     => '0',
    'population'    => 'Population',
    'paired'        => '0',
    'reads'         => '0',
    'project'       => 'HipSci_QC1_Pilot',
    'library'       => '9274735028_1559363',
    'lane_id'       => '4',
    'individual'    => 'hipscigt5466711',
    'platform'      => 'SLX',
    'center_name'   => 'SC',
    'sample'        => '271298_B03_hipscigt5466711',
    'expected_md5'  => 'd5da80818260adee189cc61acea8a9d4',
    'study'         => '2520',
    'lane'          => '9274735028_R06C01',
    'species'       => 'Homo sapiens',
    'insert_size'   => '0',
    'storage_path'  => '/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/b37c6db3-c2bb-40db-aa1b-48340681ebe1_pilot_genotyping.fcr.txt.gz'
  },
  'metadata correct for one of the genotype files';

finish;
