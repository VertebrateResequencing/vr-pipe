#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;
use Data::Dumper;

BEGIN {
    use Test::Most tests => 9;
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
foreach my $sql (VRPipe::File->create(path => file(qw(t data vrtrack_hipsci_qc1_genotyping.sql))->absolute)->slurp) {
    print $mysqlfh $sql;
}
close($mysqlfh);

# (for testing purposes in another test script, the vrtrack db here has 8
#  samples, 4 per individual, but the library name for lane 9300870057_R02C02
#  (sample fpdr_3, individual 6d3d2acf-29a5-41a2-8992-1414706a527d) was changed
#  from 283163_H01_qc1hip5533833 to 283163_G01_qc1hip5529689 (for an unrelated
#  individual))

# setup pipeline
my $output_dir = get_output_dir('genome_studio_import_genotype_files');
my $irods_dir = dir($output_dir, 'irods_import')->stringify;

# setup vrtrack datasource
ok my $ds = VRPipe::DataSource->create(
    type    => 'vrtrack',
    method  => 'analysis_genome_studio',
    source  => $ENV{VRPIPE_VRTRACK_TESTDB},
    options => { local_root_dir => $irods_dir }
  ),
  'could create a vrtrack datasource';

# check correct number of gtc file retrieved
my $results = 0;
foreach my $element (@{ get_elements($ds) }) {
    $results++;
}
is $results, 8, 'got correct number of gtc files from the vrtrack db';

# check pipeline has correct steps
ok my $pipeline = VRPipe::Pipeline->create(name => 'genome_studio_import_genotype_files'), 'able to get the genome_studio_import_genotype_files pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(irods_get_files_by_basename split_genome_studio_genotype_files genome_studio_fcr_to_vcf)], 'the pipeline has the correct steps';

my $external_gzip_file        = file(qw(t data hipsci_genotyping.fcr.txt.gz))->absolute->stringify;
my $external_reheader_penncnv = file(qw(t data reheader_penncnv.txt))->absolute->stringify;
my $snp_manifest_file         = file(qw(t data hipsci_genotyping.snp.manifest))->absolute->stringify;

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
        snp_manifest       => $snp_manifest_file,
        cleanup            => 0
    }
);

# get arrays of output files
my @irods_files;
my @lanes = qw(name 283163_A01_qc1hip5529683);
foreach my $lane (@lanes) {
    push(@irods_files, file($irods_dir, $lane . '.gtc'));
}

my @genotype_files;
my @vcf_files;
my $element_id = 0;
foreach my $sample (qw(fpdj fpdj_2 fpdj_1 fpdr fpdr_2 fpdr_1 fpdj_3 fpdr_3)) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@genotype_files, file(@output_subdirs, '2_split_genome_studio_genotype_files', $sample . '.genotyping.fcr.txt'));
    push(@vcf_files,      file(@output_subdirs, '3_genome_studio_fcr_to_vcf',           $sample . '.genotyping.fcr.vcf'));
}

#run pipeline and check outputs
ok handle_pipeline(@genotype_files, @vcf_files), 'bam import from irods and split genome studio genotype file pipeline ran ok';

#check genotype file metadata
my $meta = VRPipe::File->get(path => $genotype_files[0])->metadata;
is_deeply $meta,
  {
    'analysis_uuid' => '12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',
    'bases'         => '0',
    'withdrawn'     => '0',
    'population'    => 'Population',
    'paired'        => '0',
    'reads'         => '0',
    'project'       => 'G0325 [coreex] Wellcome Trust Strategic Award application - HIPS',
    'library'       => '283163_D02_qc1hip5533824',
    'library_tag'   => '9300870057_R06C02',
    'lane_id'       => '151',
    'individual'    => 'ca04b23b-c5b0-4389-95a3-5c7c8e6d51f2',
    'platform'      => 'SLX',
    'center_name'   => 'SC',
    'sample'        => 'fpdj',
    'expected_md5'  => 'a965097662c37a2d74a7ca429f58f726',
    'study'         => '2624',
    'control'       => 'Control',
    'lane'          => '9300870057_R06C02',
    'species'       => 'Homo sapiens',
    'insert_size'   => '0',
    'storage_path'  => '/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz'
  },
  'metadata correct for one of the genotype files';

# check the VCF is correct
is_deeply [$vcf_files[0]->slurp], [file(qw(t data hipsci_genotyping.fpdj.snp.vcf))->slurp], 'VCF file produced was as expected';
$meta = VRPipe::File->get(path => $vcf_files[0])->metadata;
is $meta->{individual}, 'ca04b23b-c5b0-4389-95a3-5c7c8e6d51f2', 'the VCF file has individual metadata';

finish;
