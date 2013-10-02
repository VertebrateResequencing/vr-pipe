#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;
use Data::Dumper;

BEGIN {
    use Test::Most tests => 11;
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
is $results, 1, 'got correct number of gtc files from the vrtrack db';

#check pipeline has correct steps
ok my $pipeline = VRPipe::Pipeline->create(name => 'genome_studio_import_genotype_files'), 'able to get the genome_studio_import_genotype_files pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(irods_get_files_by_basename split_genome_studio_genotype_files)], 'the genotype import pipeline has the correct steps';

#create external genotype gzip file for testing to override the path in gtc file metadata
my $external_gzip_source = file(qw(t data hipsci_genotyping.fcr.txt.gz));
my $gzip_dir = dir($output_dir, 'external_gzip');
$pipeline->make_path($gzip_dir);
my $external_gzip_file = file($gzip_dir, 'hipsci_genotyping.fcr.txt.gz')->stringify;
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
my @lanes = qw(name 283163_A01_qc1hip5529683);
foreach my $lane (@lanes) {
    push(@irods_files, file($irods_dir, $lane . '.gtc'));
}

my @genotype_files;
my $element_id = 0;
foreach my $sample (qw(FS18.A)) {
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
    'population'    => 'Population',
    'withdrawn'     => '0',
    'bases'         => '0',
    'analysis_uuid' => '12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',
    'paired'        => '0',
    'library_tag'   => 'R01C01',
    'reads'         => '0',
    'project'       => 'Wellcome Trust Strategic Award application – HIPS',
    'library'       => '283163_A01_qc1hip5529683',
    'lane_id'       => '58',
    'individual'    => '6d3d2acf-29a5-41a2-8992-1414706a527d',
    'sample'        => 'FS18.A',
    'center_name'   => 'SC',
    'platform'      => 'SLX',
    'study'         => '2624',
    'expected_md5'  => 'd7e10a49be4e8b1e42fe71bc68e93856',
    'lane'          => '9300870057_R01C01',
    'control'       => 'Stem cell',
    'species'       => 'Homo sapiens',
    'insert_size'   => '0',
    'storage_path'  => '/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz'
  },
  'metadata correct for one of the genotype files';

#Run penncnv pipeline using the output genotype files from the import:
$output_dir = get_output_dir('penncnv_cnv_calling');

#check pipeline has correct steps
ok my $penn_pipeline = VRPipe::Pipeline->create(name => 'penncnv_cnv_calling'), 'able to get the penncnv_cnv_calling pipeline';
my @sp_names;
foreach my $stepmember ($penn_pipeline->step_members) {
    push(@sp_names, $stepmember->step->name);
}
is_deeply \@sp_names, [qw(penncnv_detect_cnv penncnv_filter_cnv)], 'the penncnv_cnv_calling pipeline has the correct steps';

my $detect_cnv_script = '/software/vertres/bin-external/PennCNV/detect_cnv.pl';
my $detect_cnv_hmm    = '/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/genotyping/PennCNV/lib/custom.hmm';
my $detect_cnv_pfb    = '/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/genotyping/PennCNV/lib/HumanExome12v1.1.hg19.pfb';
my $filter_cnv_script = '/software/vertres/bin-external/PennCNV/filter_cnv.pl';
my $filter_numsnps    = 8;
my $filter_length     = '120k';
my $filter_confidence = 8;

VRPipe::PipelineSetup->create(
    name       => 'penncnv_calling',
    pipeline   => $penn_pipeline,
    datasource => VRPipe::DataSource->create(
        type   => 'vrpipe',
        method => 'all',
        source => 'gtc import and qc[2]',
    ),
    output_root => $output_dir,
    options     => {
        detect_cnv_script => $detect_cnv_script,
        detect_cnv_hmm    => $detect_cnv_hmm,
        detect_cnv_pfb    => $detect_cnv_pfb,
        filter_cnv_script => $filter_cnv_script,
        filter_numsnps    => $filter_numsnps,
        filter_length     => $filter_length,
        filter_confidence => $filter_confidence,
    }
);

#Get array of output files and check outputs as the pipeline is run
my @penncnv_files;
foreach my $sample (qw(FS18.A)) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 2);
    push(@penncnv_files, file(@output_subdirs, '2_penncnv_filter_cnv', $sample . '.genotyping.fcr.txt.rawcnv.filtercnv'));
}
ok handle_pipeline(@penncnv_files), 'penncnv_cnv_calling pipeline ran ok and produced the expected output files';

#check cnv file metadata
my $penn_meta = VRPipe::File->get(path => $penncnv_files[0])->metadata;
is_deeply $penn_meta,
  {
    'population'        => 'Population',
    'withdrawn'         => '0',
    'bases'             => '0',
    'analysis_uuid'     => '12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',
    'paired'            => '0',
    'library_tag'       => 'R01C01',
    'reads'             => '0',
    'project'           => 'Wellcome Trust Strategic Award application – HIPS',
    'library'           => '283163_A01_qc1hip5529683',
    'lane_id'           => '58',
    'individual'        => '6d3d2acf-29a5-41a2-8992-1414706a527d',
    'sample'            => 'FS18.A',
    'center_name'       => 'SC',
    'platform'          => 'SLX',
    'study'             => '2624',
    'expected_md5'      => 'd7e10a49be4e8b1e42fe71bc68e93856',
    'lane'              => '9300870057_R01C01',
    'control'           => 'Stem cell',
    'species'           => 'Homo sapiens',
    'insert_size'       => '0',
    'cnv_analysis_type' => 'penncnv',
    'storage_path'      => '/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz'
  },
  'metadata correct for one of the penncnv files';

finish;
