#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;
use Data::Dumper;

BEGIN {
    use Test::Most tests => 8;
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
my @sql = VRPipe::File->create(path => file(qw(t data vrtrack_hipsci_qc1_genotyping.sql))->absolute)->slurp;
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
    method  => 'analysis_genome_studio',
    source  => $ENV{VRPIPE_VRTRACK_TESTDB},
    options => { local_root_dir => $irods_dir }
  ),
  'could create a vrtrack datasource';

#check correct number of gtc file retrieved
my $results = 0;
foreach my $element (@{ get_elements($ds) }) {
    $results++;
}
is $results, 20, 'got correct number of gtc files from the vrtrack db';

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

#create external reheader file for penncnv analyses
my $reheader_penncnv = file(qw(t data reheader_penncnv.txt));
my $reheader_dir = dir($output_dir, 'reheader');
$pipeline->make_path($reheader_dir);
my $external_reheader_penncnv = file($reheader_dir, 'penncnv_reheader.txt')->stringify;
copy($reheader_penncnv, $external_reheader_penncnv);

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
        reheader_penncnv   => $external_reheader_penncnv,        
        cleanup            => 1
    }
);

#get arrays of output files
my @irods_files;
my @lanes = qw(9300870166_R06C01 9300870166_R05C01 9300870166_R04C01 9300870166_R03C01 9300870166_R02C02 9300870166_R02C01 9300870166_R01C02 9300870166_R01C01 9300870057_R06C02 9300870057_R06C01 9300870057_R05C02 9300870057_R05C01 9300870057_R04C02 9300870057_R04C01 9300870057_R03C02 9300870057_R03C01 9300870057_R02C02 9300870057_R02C01 9300870057_R01C02 9300870057_R01C01);
foreach my $lane (@lanes) {
    push(@irods_files, file($irods_dir, $lane . '.gtc'));
}

my @genotype_files;
my $element_id = 0;
foreach my $sample (qw(qc1hip5533830 qc1hip5533829 qc1hip5533828 qc1hip5533827 qc1hip5533832 qc1hip5533826 qc1hip5533831 qc1hip5533825 qc1hip5533824 qc1hip5529688 qc1hip5533823 qc1hip5529687 qc1hip5533822 qc1hip5529686 qc1hip5533821 qc1hip5529685 qc1hip5533833 qc1hip5529684 qc1hip5529689 qc1hip5529683)) {
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
    'analysis_uuid' => '12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',
    'bases' => '0',
    'withdrawn' => '0',
    'population' => 'Population',
    'paired' => '0',
    'reads' => '0',
    'project' => 'Wellcome Trust Strategic Award application – HIPS',
    'library' => '283163_B03_qc1hip5533830',
    'lane_id' => '1',
    'individual' => '2a39941c-12b2-41bf-92f3-70b88b66a3a4',
    'platform' => 'SLX',
    'center_name' => 'SC',
    'sample' => 'qc1hip5533830',
    'expected_md5' => '7793f115dadaa5e0a2b4aae5aca89ce9',
    'study' => '2624',
    'lane' => '9300870166_R06C01',
    'species' => 'Homo sapiens',
    'insert_size' => '0',
    'storage_path' => '/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz'
  },
  'metadata correct for one of the genotype files';

#Run penncnv pipeline using the output genotype files from the import:
#Add test code here!
$output_dir = get_output_dir('penncnv_analysis');
VRPipe::PipelineSetup->create(
    name       => 'penncnv_calling',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'gtc import and qc[2]',

    ),
    output_root => $output_dir,
    pipeline    => VRPipe::Pipeline->create(name => 'penncnv'), #whatever name of pipeline is.....
    options     => {
		#options go here.....
    }
);

ok handle_pipeline(), 'penncnv pipeline ran';

finish;
