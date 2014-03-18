#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 9;
    use VRPipeTest;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_AUTHOR_TESTS WAREHOUSE_DATABASE WAREHOUSE_HOST WAREHOUSE_PORT WAREHOUSE_USER)],
        required_exe => [qw(imeta)]
    );
    use TestPipelines;
}

# these are author-only tests (Sanger-specific) to make sure that we get all
# the desired metadata out of irods and warehouse for all the different types
# of data we process. We'll also test that we can then take this VRPipe metadata
# on files and store all the tracking information in an external database

my $output_root = get_output_dir('datasource_irods_import_dir');

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

my $files = 0;
foreach my $result (map { result_with_inflated_paths($_) } @{ get_elements($ds) }) {
    $files++;
    next if $files > 1;
    my ($path) = @{ $result->{paths} };
    my $file = VRPipe::File->create(path => $path);
    is_deeply $file->metadata,
      {
        'analysis_uuid'        => 'cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',
        'sample_cohort'        => '6d3d2acf-29a5-41a2-8992-1414706a527d',
        'study_id'             => [2622, 2625],
        'study_title'          => ['G0325 [collection qc1] Wellcome Trust Strategic Award application – HIPS', 'G0325 [gex] Wellcome Trust Strategic Award application – HIPS'],
        'sample'               => 'qc1hip5529781',
        'public_name'          => 'fpdr',
        'sample_created_date'  => '2013-05-10 06:45:32',
        'beadchip_section'     => 'E',
        'sample_consent'       => '1',
        'taxon_id'             => '9606',
        'irods_path'           => '/archive/GAPI/exp/infinium/41/b6/f3/9252616016_E_Grn.idat',
        'sample_common_name'   => 'Homo Sapien',
        'sample_control'       => '1',
        'beadchip'             => '9252616016',
        'sample_id'            => '1625281',
        'md5'                  => '41b6f345a2f92f095f09a7ff22bbbc00',
        'sample_supplier_name' => 'face6d88-7e90-4215-aa80-fb2c3df5a4ed'
      },
      'correct file metadata was present on the first idat file';
}
is $files, 7, 'idat datasource returned the correct number of files';

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

$files = 0;
foreach my $result (map { result_with_inflated_paths($_) } @{ get_elements($ds) }) {
    $files++;
    next if $files > 1;
    my ($path) = @{ $result->{paths} };
    my $file = VRPipe::File->create(path => $path);
    is_deeply $file->metadata,
      {
        'analysis_uuid'        => ['45a53a77-50bc-4062-b9cb-8dfe82e589f2', '12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab', '7da360df-218d-4cd0-994d-37fefd0ffaa9', '3f5acca0-304c-480f-8a61-3e68c33c707d'],
        'infinium_well'        => 'F01',
        'study_id'             => '2624',
        'sample_cohort'        => '27af9a9b-01b2-4cb6-acef-ea52d83e3d26',
        'infinium_plate'       => 'WG0206884-DNA',
        'infinium_sample'      => '283163_F01_qc1hip5529688',
        'public_name'          => 'fpdk_3',
        'sample'               => 'qc1hip5529688',
        'beadchip_design'      => 'HumanCoreExome-12v1-0',
        'sample_created_date'  => '2013-05-10 06:44:46',
        'sample_consent'       => '1',
        'beadchip_section'     => 'R06C01',
        'irods_path'           => '/archive/GAPI/gen/infinium/17/b7/15/9300870057_R06C01.gtc',
        'taxon_id'             => '9606',
        'sample_control'       => '0',
        'sample_common_name'   => 'Homo Sapien',
        'study_title'          => 'G0325 [coreex] Wellcome Trust Strategic Award application – HIPS',
        'beadchip'             => '9300870057',
        'sample_id'            => '1625188',
        'md5'                  => '17b7159554bca4ff4376384b385da51f',
        'sample_supplier_name' => '87e7ee6f-e16f-41f6-94c5-194933e2b192'
      },
      'correct file metadata was present on the first gtc file';
}
is $files, 20, 'gtc datasource returned the correct number of files';

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

exit;
