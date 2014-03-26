#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 21;
    use VRPipeTest;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB VRPIPE_AUTHOR_TESTS WAREHOUSE_DATABASE WAREHOUSE_HOST WAREHOUSE_PORT WAREHOUSE_USER)],
        required_exe => [qw(imeta)]
    );
    use TestPipelines;
    
    use_ok('VRTrack::Factory');
}

# these are author-only tests (Sanger-specific) to make sure that we get all
# the desired metadata out of irods and warehouse for all the different types
# of data we process. We'll also test that we can then take this VRPipe metadata
# on files and store all the tracking information in an external database

my $output_root = get_output_dir('irods_metadata_populator');

# idat expression files
ok my $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'archive',
    options => {
        file_query     => q[study_id = 2625 and type = idat and dcterms:created '<' 2013-06-07],
        local_root_dir => $output_root
    }
  ),
  'could create an irods datasource for idat files';

create_fresh_vrtrack_test_db();
ok my $pipeline = VRPipe::Pipeline->create(name => 'vrtrack_populate_from_irods_and_download_files'), 'able to get the vrtrack_populate_from_irods_and_download_files pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(vrtrack_populate_from_vrpipe_metadata irods_analysis_files_download)], 'the pipeline has the correct steps';

VRPipe::PipelineSetup->create(
    name        => 'idat populate',
    datasource  => $ds,
    output_root => $output_root,
    pipeline    => $pipeline,
    options     => {
        vrtrack_db         => $ENV{VRPIPE_VRTRACK_TESTDB},
        vrlane_storage_dir => $output_root
    }
);

my @expected_output_files = (file($output_root, '/archive/GAPI/exp/analysis/69/61/88/hipsci_7samples_2013-05-15/hipsci_2013-05-15_Sample_Probe_Profile.txt.gz'), file($output_root, '/archive/GAPI/exp/analysis/69/61/88/hipsci_7samples_2013-05-15/hipsci_2013-05-15_annotation.txt'));
ok handle_pipeline(@expected_output_files), 'vrtrack_populate_from_vrpipe_metadata pipeline ran ok for idat files';

my $files = 0;
foreach my $result (map { result_with_inflated_paths($_) } @{ get_elements($ds) }) {
    $files++;
    next if $files > 1;
    my ($path) = @{ $result->{paths} };
    my $file = VRPipe::File->create(path => $path);
    is_deeply $file->metadata,
      {
        'analysis_uuid'           => 'cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',
        'sample_cohort'           => '6d3d2acf-29a5-41a2-8992-1414706a527d',
        'study_id'                => 2625,
        'study_title'             => 'G0325 [gex] Wellcome Trust Strategic Award application – HIPS',
        'sample'                  => 'qc1hip5529781',
        'public_name'             => 'fpdr',
        'sample_created_date'     => '2013-05-10 06:45:32',
        'beadchip_section'        => 'E',
        'sample_consent'          => '1',
        'taxon_id'                => '9606',
        'irods_path'              => '/archive/GAPI/exp/infinium/41/b6/f3/9252616016_E_Grn.idat',
        'sample_common_name'      => 'Homo Sapien',
        'sample_control'          => '1',
        'beadchip'                => '9252616016',
        'sample_id'               => '1625281',
        'md5'                     => '41b6f345a2f92f095f09a7ff22bbbc00',
        'sample_supplier_name'    => 'face6d88-7e90-4215-aa80-fb2c3df5a4ed',
        'irods_analysis_files'    => ['/archive/GAPI/exp/analysis/69/61/88/hipsci_7samples_2013-05-15/hipsci_2013-05-15_Sample_Probe_Profile.txt.gz', '/archive/GAPI/exp/analysis/69/61/88/hipsci_7samples_2013-05-15/hipsci_2013-05-15_annotation.txt',],
        'irods_local_storage_dir' => $output_root
      },
      'correct file metadata was present on the first idat file';
}
is $files, 7, 'idat datasource returned the correct number of files';

my $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
my @lanes = $vrtrack->get_lanes();
is @lanes, 7, 'the correct number of idat lanes were populated';
my $lane_info = lane_info(grep { $_->name eq '9252616016_E_Grn' } @lanes);

my $expected = {
    lane_name             => '9252616016_E_Grn',
    raw_reads             => undef,
    is_paired             => undef,
    lane_acc              => 'cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',
    storage_path          => file($output_root, 'archive/GAPI/exp/analysis/69/61/88/hipsci_7samples_2013-05-15')->stringify,
    file_name             => '9252616016_E_Grn.idat',
    file_type             => 8,
    file_md5              => '41b6f345a2f92f095f09a7ff22bbbc00',
    library_name          => 'qc1hip5529781.9252616016.E',
    library_ssid          => 6016281,
    library_ts            => '9252616016_E',
    project_id            => 2625,
    project_name          => 'G0325 [gex] Wellcome Trust Strategic Award application - HIPS',
    study_acc             => 2625,
    sample_name           => 'fpdr',
    sample_hierarchy_name => 'face6d88-7e90-4215-aa80-fb2c3df5a4ed',
    sample_ssid           => '1625281',
    species_taxon_id      => 9606,
    species_name          => 'Homo Sapien',
    individual_name       => '6d3d2acf-29a5-41a2-8992-1414706a527d',
    individual_acc        => undef
};

is_deeply $lane_info, $expected, 'VRTrack was correctly populated for the first idat lane';

# gtc genotyping files
ok $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'archive',
    options => {
        file_query     => q[study_id = 2624 and type = gtc and dcterms:created '<' 2013-06-01],
        local_root_dir => $output_root
    }
  ),
  'could create an irods datasource for gtc files';

create_fresh_vrtrack_test_db();
VRPipe::PipelineSetup->create(
    name        => 'gtc populate',
    datasource  => $ds,
    output_root => $output_root,
    pipeline    => $pipeline,
    options     => {
        vrtrack_db         => $ENV{VRPIPE_VRTRACK_TESTDB},
        vrlane_storage_dir => $output_root,
    }
);

@expected_output_files = (file($output_root, '/archive/GAPI/gen/analysis/be/f4/37/coreex_hips/20130808/coreex_hips_20130808.fcr.txt.gz'));
ok handle_pipeline(@expected_output_files), 'vrtrack_populate_from_vrpipe_metadata pipeline ran ok for gtc files';

$files = 0;
foreach my $result (map { result_with_inflated_paths($_) } @{ get_elements($ds) }) {
    $files++;
    next if $files > 1;
    my ($path) = @{ $result->{paths} };
    my $file = VRPipe::File->create(path => $path);
    is_deeply $file->metadata,
      {
        'analysis_uuid'           => ['45a53a77-50bc-4062-b9cb-8dfe82e589f2', '12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab', '7da360df-218d-4cd0-994d-37fefd0ffaa9', '3f5acca0-304c-480f-8a61-3e68c33c707d'],
        'infinium_well'           => 'F01',
        'study_id'                => '2624',
        'sample_cohort'           => '27af9a9b-01b2-4cb6-acef-ea52d83e3d26',
        'infinium_plate'          => 'WG0206884-DNA',
        'infinium_sample'         => '283163_F01_qc1hip5529688',
        'public_name'             => 'fpdk_3',
        'sample'                  => 'qc1hip5529688',
        'beadchip_design'         => 'HumanCoreExome-12v1-0',
        'sample_created_date'     => '2013-05-10 06:44:46',
        'sample_consent'          => '1',
        'beadchip_section'        => 'R06C01',
        'irods_path'              => '/archive/GAPI/gen/infinium/17/b7/15/9300870057_R06C01.gtc',
        'taxon_id'                => '9606',
        'sample_control'          => '0',
        'sample_common_name'      => 'Homo Sapien',
        'study_title'             => 'G0325 [coreex] Wellcome Trust Strategic Award application – HIPS',
        'beadchip'                => '9300870057',
        'sample_id'               => '1625188',
        'md5'                     => '17b7159554bca4ff4376384b385da51f',
        'sample_supplier_name'    => '87e7ee6f-e16f-41f6-94c5-194933e2b192',
        'irods_analysis_files'    => '/archive/GAPI/gen/analysis/be/f4/37/coreex_hips/20130808/coreex_hips_20130808.fcr.txt.gz',
        'irods_local_storage_dir' => $output_root
      },
      'correct file metadata was present on the first gtc file';
}
is $files, 20, 'gtc datasource returned the correct number of files';

$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
@lanes = $vrtrack->get_lanes();
is @lanes, 20, 'the correct number of gtc lanes were populated';
$lane_info = lane_info(grep { $_->name eq '9300870057_R06C01' } @lanes);

$expected = {
    lane_name             => '9300870057_R06C01',
    raw_reads             => undef,
    is_paired             => undef,
    lane_acc              => '3f5acca0-304c-480f-8a61-3e68c33c707d',
    storage_path          => file($output_root, 'archive/GAPI/gen/analysis/be/f4/37/coreex_hips/20130808/coreex_hips_20130808.fcr.txt.gz')->stringify,
    file_name             => '9300870057_R06C01.gtc',
    file_type             => 7,
    file_md5              => '17b7159554bca4ff4376384b385da51f',
    library_name          => '283163_F01_qc1hip5529688',
    library_ssid          => 57188,
    library_ts            => '9300870057_R06C01',
    project_id            => 2624,
    project_name          => 'G0325 [coreex] Wellcome Trust Strategic Award application - HIPS',
    study_acc             => 2624,
    sample_name           => 'fpdk_3',
    sample_hierarchy_name => '87e7ee6f-e16f-41f6-94c5-194933e2b192',
    sample_ssid           => '1625188',
    species_taxon_id      => 9606,
    species_name          => 'Homo Sapien',
    individual_name       => '27af9a9b-01b2-4cb6-acef-ea52d83e3d26',
    individual_acc        => undef
};

is_deeply $lane_info, $expected, 'VRTrack was correctly populated for the first gtc lane';

# bam sequencing data
ok $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'seq',
    options => {
        file_query     => q[study_id = 2547 and type = bam | grep -v phix | grep -v '#0'],
        local_root_dir => $output_root
    }
  ),
  'could create an irods datasource for bam files';

create_fresh_vrtrack_test_db();
VRPipe::PipelineSetup->create(
    name        => 'bam populate',
    datasource  => $ds,
    output_root => $output_root,
    pipeline    => $pipeline,
    options     => {
        vrtrack_db                 => $ENV{VRPIPE_VRTRACK_TESTDB},
        vrlane_storage_dir         => $output_root,
        irods_download_input_files => 1
    }
);

@expected_output_files = (file($output_root, '/seq/9417/9417_4#1.bam'), file($output_root, '/seq/9417/9417_4#2.bam'), file($output_root, '/seq/9417/9417_4#3.bam'), file($output_root, '/seq/9417/9417_4#4.bam'));
ok handle_pipeline(@expected_output_files), 'vrtrack_populate_from_vrpipe_metadata pipeline ran ok for bam files';

$files = 0;
foreach my $result (map { result_with_inflated_paths($_) } @{ get_elements($ds) }) {
    $files++;
    next if $files > 1;
    my ($path) = @{ $result->{paths} };
    my $file = VRPipe::File->create(path => $path);
    is_deeply $file->metadata,
      {
        'study_id'                => '2547',
        'is_paired_read'          => '1',
        'library_id'              => '6784051',
        'ebi_sub_acc'             => 'ERA214806',
        'library'                 => 'MEK_res_1 6784051',
        'study_title'             => 'De novo and acquired resistance to MEK inhibitors ',
        'target'                  => '1',
        'reference'               => '/lustre/scratch109/srpipe/references/Mus_musculus/GRCm38/all/bwa/Mus_musculus.GRCm38.68.dna.toplevel.fa',
        'alignment'               => '1',
        'sample'                  => 'MEK_res_1',
        'sample_accession_number' => 'ERS215816',
        'tag'                     => 'AACGTGAT',
        'study'                   => 'De novo and acquired resistance to MEK inhibitors',
        'lane'                    => '4',
        'sample_created_date'     => '2013-02-21 07:51:47',
        'ebi_sub_md5'             => '675b0b2b2f5991aa5a4695bb1914c0c7',
        'ebi_run_acc'             => 'ERR274701',
        'taxon_id'                => '10090',
        'irods_path'              => '/seq/9417/9417_4#1.bam',
        'ebi_sub_date'            => '2013-05-26',
        'total_reads'             => '70486376',
        'id_run'                  => '9417',
        'tag_index'               => '1',
        'sample_common_name'      => 'Mouse',
        'study_accession_number'  => 'ERP002262',
        'manual_qc'               => '1',
        'sample_id'               => '1571700',
        'md5'                     => '675b0b2b2f5991aa5a4695bb1914c0c7'
      },
      'correct file metadata was present on the first bam file';
}
is $files, 4, 'bam datasource returned the correct number of files';

$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
@lanes = $vrtrack->get_lanes();
is @lanes, 4, 'the correct number of bam lanes were populated';
$lane_info = lane_info(grep { $_->name eq '9417_4#1' } @lanes);

$expected = {
    lane_name             => '9417_4#1',
    raw_reads             => 70486376,
    is_paired             => 1,
    lane_acc              => undef,
    storage_path          => undef,
    file_name             => '9417_4#1.bam',
    file_type             => 4,
    file_md5              => '675b0b2b2f5991aa5a4695bb1914c0c7',
    library_name          => 'MEK_res_1 6784051',
    library_ssid          => 6784051,
    library_ts            => undef,
    project_id            => 2547,
    project_name          => 'De novo and acquired resistance to MEK inhibitors ',
    study_acc             => 'ERP002262',
    sample_name           => 'MEK_res_1',
    sample_hierarchy_name => 'MEK_res_1',
    sample_ssid           => '1571700',
    species_taxon_id      => 10090,
    species_name          => 'Mouse',
    individual_name       => 'MEK_res_1',
    individual_acc        => 'ERS215816'
};

is_deeply $lane_info, $expected, 'VRTrack was correctly populated for the first bam lane';

exit;

sub create_fresh_vrtrack_test_db {
    my %cd = VRTrack::Factory->connection_details('rw');
    open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
    print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    foreach my $sql (VRTrack::VRTrack->schema()) {
        print $mysqlfh $sql;
    }
    close($mysqlfh);
}

sub lane_info {
    my $lane   = shift;
    my ($file) = @{ $lane->files };
    my %h      = $vrtrack->lane_hierarchy_objects($lane);
    
    my $info = {
        lane_name             => $lane->name,
        raw_reads             => $lane->raw_reads,
        is_paired             => $lane->is_paired,
        lane_acc              => $lane->acc,
        storage_path          => $lane->storage_path,
        file_name             => $file->name,
        file_type             => $file->type,
        file_md5              => $file->md5,
        library_name          => $h{library}->name,
        library_ssid          => $h{library}->ssid,
        library_ts            => $h{library}->library_tag_sequence,
        project_id            => $h{project}->ssid,
        project_name          => $h{project}->name,
        study_acc             => $h{study}->acc,
        sample_name           => $h{sample}->name,
        sample_hierarchy_name => $h{sample}->hierarchy_name,
        sample_ssid           => $h{sample}->ssid,
        species_taxon_id      => $h{species}->taxon_id,
        species_name          => $h{species}->name,
        individual_name       => $h{individual}->name,
        individual_acc        => $h{individual}->acc
    };
    
    return $info;
}
