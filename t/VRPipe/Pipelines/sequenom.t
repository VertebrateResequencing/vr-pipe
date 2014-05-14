#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 20;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_AUTHOR_TESTS)],
    );
    use TestPipelines;
}

ok my $pipeline = VRPipe::Pipeline->create(name => 'sequenom_import_from_irods_and_covert_to_vcf'), 'able to get the sequenom_import_from_irods_and_covert_to_vcf pipeline';
my @step_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@step_names, $stepmember->step->name);
}
is_deeply \@step_names, [qw(irods_get_files_by_basename sequenom_csv_to_vcf)], 'the sequenom_import_from_irods_and_covert_to_vcf pipeline has the correct steps';

my $output_root = get_output_dir('sequenom_import');
my $sequenom_plex_storage_dir = dir($output_root, 'sequenom_plex_storage_dir');
mkdir($sequenom_plex_storage_dir);

my $file_query = q[sequenom_plate LIKE '%' and study_id = 2622 and dcterms:created '<' 2013-07-26];

ok my $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'seq',
    options => {
        file_query     => $file_query,
        local_root_dir => $output_root
    }
  ),
  'could create the irods datasource';

my $ps = VRPipe::PipelineSetup->create(
    name        => 'sequenom import',
    datasource  => $ds,
    output_root => $output_root,
    pipeline    => $pipeline,
    options     => { sequenom_plex_storage_dir => $sequenom_plex_storage_dir }
);

my (@vcf_files, @index_files);
my $element_id = 1;
foreach my $basname (qw(QC288261____20130701_G01 QC288261____20130701_C01 QC288261____20130701_A01 QC288261____20130701_E01)) {
    my @output_subdirs = output_subdirs($element_id++, 1);
    push(@vcf_files, file(@output_subdirs, '2_sequenom_csv_to_vcf', $basname . '.vcf.gz'));
    push(@index_files, $vcf_files[-1] . '.csi');
}
ok handle_pipeline(@vcf_files, @index_files), 'sequenom_import_from_irods_and_covert_to_vcf pipeline ran ok and produced the expected output files';

# this is already tested in another script, but we'll do an end-to-end test with
# the real-world data to ensure we can confirm the genotypes are sensible
$output_root = get_output_dir('vcf_merge_and_compare_genotypes');
$ds          = VRPipe::DataSource->create(
    type    => 'vrpipe',
    method  => 'group_by_metadata',
    source  => '1[2:vcf_files]',
    options => { metadata_keys => 'sample_cohort' }
);

$pipeline = VRPipe::Pipeline->create(name => 'vcf_merge_and_compare_genotypes');

$ps = VRPipe::PipelineSetup->create(
    name        => 'genotype comparison',
    datasource  => $ds,
    output_root => $output_root,
    pipeline    => $pipeline,
    options     => {}
);

my (@merged_vcf_files, @gtypex_files, @expected_metadata);
foreach my $element_id (5, 6, 7) {
    my @output_subdirs = output_subdirs($element_id, 2);
    push(@merged_vcf_files, file(@output_subdirs, '1_vcf_merge_different_samples', 'merged.vcf.gz'));
    push(@gtypex_files,     file(@output_subdirs, '2_vcf_genotype_comparison',     'merged.vcf.gz.gtypex'));
    
    my $individual = VRPipe::DataElement->get(id => $element_id)->metadata->{group};
    my %expected = (sample_cohort => $individual);
    if ($individual eq '20f8a331-69ac-4510-94ab-e3a69c50e46f') {
        $expected{genotype_maximum_deviation} = "0.000000:HPSI0813i-ffdb_4_QC1Hip-2";
        $expected{sequenom_gender}            = 'M';
        $expected{sample}                     = undef;
        $expected{public_name}                = undef;
    }
    elsif ($individual eq '3d52354f-8d84-457d-a668-099a758f0e7b') {
        $expected{genotype_maximum_deviation} = '0.000000:lofv_33_QC1Hip-4';
        $expected{sequenom_gender}            = 'F';
        $expected{sample}                     = 'QC1Hip-4';
        $expected{public_name}                = 'lofv_33';
    }
    else {
        $expected{genotype_maximum_deviation} = "0.000000:HPSI0813i-ffdc_5_QC1Hip-3";
        $expected{sequenom_gender}            = 'M';
        $expected{sample}                     = 'QC1Hip-3';
        $expected{public_name}                = 'HPSI0813i-ffdc_5';
    }
    push(@expected_metadata, \%expected);
}

ok handle_pipeline(@merged_vcf_files, @gtypex_files), 'vcf_merge_and_compare_genotypes pipeline ran ok';

foreach my $vcf_path (@merged_vcf_files) {
    my $meta = VRPipe::File->get(path => $vcf_path)->metadata;
    my $expected = shift @expected_metadata;
    foreach my $key (qw(sample_cohort genotype_maximum_deviation sequenom_gender sample public_name)) {
        is $meta->{$key}, $expected->{$key}, "$key metadata was correct for one of the merged VCF files";
    }
}

finish;
exit;
