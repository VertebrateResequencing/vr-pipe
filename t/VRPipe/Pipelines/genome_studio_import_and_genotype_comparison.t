#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;
use Data::Dumper;

BEGIN {
    use Test::Most tests => 37;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB VRPIPE_AUTHOR_TESTS WAREHOUSE_DATABASE WAREHOUSE_HOST WAREHOUSE_PORT WAREHOUSE_USER)],
        required_exe => [qw(iget iquest fcr-to-vcf sort bgzip)]
    );
    use TestPipelines;
    
    use_ok('VRTrack::Factory');
}

# create an empty vrtrack db
my %cd = VRTrack::Factory->connection_details('rw');
open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
foreach my $sql (VRTrack::VRTrack->schema()) {
    print $mysqlfh $sql;
}
close($mysqlfh);

my $output_dir = get_output_dir('genome_studio');
my $irods_dir = dir($output_dir, 'irods_import')->stringify;

# create a setup to import idat files from irods
my $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'archive',
    options => {
        file_query      => q[study_id = 2624 and type = gtc and dcterms:created '<' 2013-06-01],
        local_root_dir  => $irods_dir,
        update_interval => 99999999
    }
);

my $pipeline = VRPipe::Pipeline->create(name => 'vrtrack_populate_from_irods_and_download_files');

VRPipe::PipelineSetup->create(
    name        => 'gtc populate',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        vrtrack_db         => $ENV{VRPIPE_VRTRACK_TESTDB},
        vrlane_storage_dir => $irods_dir
    }
);

my @analysis_files = (file($irods_dir, '/archive/GAPI/gen/analysis/74/39/87/coreex_hips/20130613/coreex_hips_20130613.fcr.txt.gz'));
ok handle_pipeline(@analysis_files), 'vrtrack_populate_from_vrpipe_metadata pipeline ran ok and got the analysis files';

# create split->vcf setup using the output gtc files from the import
$output_dir = get_output_dir('genome_studio_split_and_convert_to_vcf');

# check pipeline has correct steps
ok my $split_convert_pipeline = VRPipe::Pipeline->create(name => 'genome_studio_split_and_convert_to_vcf'), 'able to get the genome_studio_split_and_convert_to_vcf pipeline';
my @s_names;
foreach my $stepmember ($split_convert_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(split_genome_studio_genotype_files illumina_coreexome_manifest_to_map genome_studio_fcr_to_vcf)], 'the pipeline has the correct steps';

my $external_reheader_penncnv = file(qw(t data reheader_penncnv.txt))->absolute->stringify;
my $manifest_file_soure       = file(qw(t data hipsci_genotyping.snp.manifest))->absolute->stringify;
my $manifest_file             = file($output_dir, 'hipsci_genotyping.manifest')->stringify;
copy($manifest_file_soure, $manifest_file);

# create pipeline setup
my $split_convert_setup = VRPipe::PipelineSetup->create(
    name       => 'gtc split and convert',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'gtc populate[2:input_files]',
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => $split_convert_pipeline,
    options     => {
        reheader_penncnv   => $external_reheader_penncnv,
        coreexome_manifest => $manifest_file,
        cleanup            => 0
    }
);

my @samples    = qw(fpdk_3_qc1hip5529688 fpdr_qc1hip5529685 fpdl_3_qc1hip5533827 fpdj_1_qc1hip5533821 fpdm_3_qc1hip5533831 fpdl_2_qc1hip5533826 fpdl_1_qc1hip5533825 fpdk_2_qc1hip5529687 fpdm_2_qc1hip5533830 fpdj_qc1hip5533824 fpdl_qc1hip5533828 fpdr_3_qc1hip5533833 fpdj_3_qc1hip5533823 fpdj_2_qc1hip5533822 fpdr_1_qc1hip5529683 fpdk_qc1hip5529689 fpdk_1_qc1hip5529686 fpdm_qc1hip5533832 fpdr_2_qc1hip5529684 fpdm_1_qc1hip5533829);
my $element_id = 20;
my @genotype_files;
my @vcf_files;
foreach my $sample (@samples) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 2);
    
    # for testing purposes we'll fake a sample swap between 2 of the cohorts by
    # altering sample names and cohort
    if ($sample eq 'fpdk_2_qc1hip5529687' || $sample eq 'fpdj_2_qc1hip5533822') {
        my ($input_file) = @{ VRPipe::DataElement->get(id => $element_id)->files };
        if ($sample eq 'fpdk_2_qc1hip5529687') {
            $input_file->add_metadata({ public_name => 'fpdj_2', sample => 'qc1hip5533822', sample_cohort => 'ca04b23b-c5b0-4389-95a3-5c7c8e6d51f2' });
            $sample = 'fpdj_2_qc1hip5533822';
        }
        else {
            $input_file->add_metadata({ public_name => 'fpdk_2', sample => 'qc1hip5529687', sample_cohort => '27af9a9b-01b2-4cb6-acef-ea52d83e3d26' });
            $sample = 'fpdk_2_qc1hip5529687';
        }
    }
    
    my ($sanger_name) = $sample =~ /(qc\S+)/;
    push(@genotype_files, file(@output_subdirs, '1_split_genome_studio_genotype_files', $sanger_name . '.genotyping.fcr.txt'));
    
    push(@vcf_files, file(@output_subdirs, '3_genome_studio_fcr_to_vcf', "$sample/$sample.vcf.gz"));
}

# run pipeline and check outputs
ok handle_pipeline(@genotype_files, @vcf_files), 'genome_studio_split_and_convert_to_vcf pipeline ran ok';

# check genotype file metadata
my $meta = VRPipe::File->get(path => $genotype_files[0])->metadata;
is_deeply $meta,
  {
    analysis_uuid           => [qw(45a53a77-50bc-4062-b9cb-8dfe82e589f2 12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab 3f5acca0-304c-480f-8a61-3e68c33c707d)],
    beadchip                => '9300870057',
    beadchip_design         => 'HumanCoreExome-12v1-0',
    beadchip_section        => 'R06C01',
    infinium_plate          => 'WG0206884-DNA',
    infinium_sample         => '283163_F01_qc1hip5529688',
    infinium_well           => 'F01',
    irods_analysis_files    => '/archive/GAPI/gen/analysis/74/39/87/coreex_hips/20130613/coreex_hips_20130613.fcr.txt.gz',
    irods_local_storage_dir => $irods_dir,
    irods_path              => '/archive/GAPI/gen/infinium/17/b7/15/9300870057_R06C01.gtc',
    md5                     => '17b7159554bca4ff4376384b385da51f',
    public_name             => 'fpdk_3',
    sample                  => 'qc1hip5529688',
    sample_cohort           => '27af9a9b-01b2-4cb6-acef-ea52d83e3d26',
    sample_common_name      => 'Homo Sapien',
    sample_consent          => 1,
    sample_control          => 0,
    sample_created_date     => '2013-05-10 06:44:46',
    sample_id               => 1625188,
    sample_supplier_name    => '87e7ee6f-e16f-41f6-94c5-194933e2b192',
    study_id                => 2624,
    study_title             => 'G0325 [coreex] Wellcome Trust Strategic Award application – HIPS',
    taxon_id                => 9606
  },
  'metadata correct for one of the genotype files';

# check the VCF is correct
is_deeply [vcf_lines($vcf_files[0])], [vcf_lines(file(qw(t data fpdk_3_qc1hip5529688.vcf.gz)))], 'VCF file produced was as expected';
$meta = VRPipe::File->get(path => $vcf_files[0])->metadata;
is $meta->{sample_cohort}, '27af9a9b-01b2-4cb6-acef-ea52d83e3d26', 'the VCF file has sample_cohort metadata';

# we'll take this opportunity to test the vcf_merge_and_compare_genotypes
# pipeline as well
$output_dir = get_output_dir('vcf_merge_and_compare_genotypes');
my $vrpipe_ds = VRPipe::DataSource->create(
    type    => 'vrpipe',
    method  => 'group_by_metadata',
    source  => '2[3:vcf_files]',
    options => { metadata_keys => 'sample_cohort|beadchip' }
);

# check pipeline has correct steps
ok my $gt_pipeline = VRPipe::Pipeline->create(name => 'vcf_merge_and_compare_genotypes'), 'able to get the vcf_merge_and_compare_genotypes pipeline';
@s_names = ();
foreach my $stepmember ($gt_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(vcf_index vcf_merge_different_samples vcf_index vcf_genotype_comparison)], 'the pipeline has the correct steps';

# create pipeline setup
VRPipe::PipelineSetup->create(
    name        => 'genotype comparison',
    datasource  => $vrpipe_ds,
    output_root => $output_dir,
    pipeline    => $gt_pipeline,
    options     => {}
);

my (@merged_vcf_files, @gtypex_files, @expected_metadata);
foreach my $element_id (41 .. 45) {
    my @output_subdirs = output_subdirs($element_id, 3);
    push(@genotype_files,   file(@output_subdirs, '2_vcf_merge_different_samples', 'merged.vcf.gz'));
    push(@merged_vcf_files, file(@output_subdirs, '4_vcf_genotype_comparison',     'merged.vcf.gz.gtypex'));
    
    my $group = VRPipe::DataElement->get(id => $element_id)->metadata->{group};
    my %expected = (group => $group);
    if ($group eq '2a39941c-12b2-41bf-92f3-70b88b66a3a4|9300870166') {
        $expected{genotype_maximum_deviation} = ['==', 0, 'fpdm_1_qc1hip5533829'];
        #$expected{genotype_maximum_deviation} = ['>', 10, 'fpdr_3'];
    }
    elsif ($group eq '27af9a9b-01b2-4cb6-acef-ea52d83e3d26|9300870057') {
        $expected{genotype_maximum_deviation} = ['>=', 18, 'fpdk_2_qc1hip5529687'];
    }
    elsif ($group eq 'ca04b23b-c5b0-4389-95a3-5c7c8e6d51f2|9300870057') {
        $expected{genotype_maximum_deviation} = ['>=', 18, 'fpdj_2_qc1hip5533822'];
    }
    elsif ($group eq '647d3009-5603-4b07-bf02-6161c8662f46|9300870166') {
        $expected{genotype_maximum_deviation} = ['==', 0, 'fpdl_qc1hip5533828'];
    }
    elsif ($group eq '6d3d2acf-29a5-41a2-8992-1414706a527d|9300870057') {
        $expected{genotype_maximum_deviation} = ['==', 0, 'fpdr_2_qc1hip5529684'];
    }
    push(@expected_metadata, \%expected);
}

ok handle_pipeline(@merged_vcf_files, @gtypex_files), 'vcf_merge_and_compare_genotypes pipeline ran ok';

foreach my $vcf_path (@merged_vcf_files) {
    $meta = VRPipe::File->get(path => $vcf_path)->metadata;
    my $expected = shift @expected_metadata;
    
    is "$meta->{sample_cohort}|$meta->{beadchip}", $expected->{group}, "sample_cohort and beadchip metadata was correct for one of the merged VCF files";
    
    my ($cmp, $eval, $esample) = @{ $expected->{genotype_maximum_deviation} };
    my ($aval, $asample) = split(':', $meta->{genotype_maximum_deviation});
    cmp_ok $aval, $cmp, $eval, "genotype_maximum_deviation metadata value was correct for one of the merged VCF files";
    is $asample, $esample, "genotype_maximum_deviation metadata sample was correct for one of the merged VCF files";
}

# we'll also take the opportunity to test bcftools cnv caller pipeline, since that
# also uses files from the genome studio import
SKIP: {
    my $num_tests = 5;
    skip "hipsci bcftools cnv calling tests disabled without bcftools_cnv_caller in your path", $num_tests if 1; #unless can_execute('bcftools');
    
    $output_dir = get_output_dir('bcftools_cnv_caller');
    
    # check pipeline has correct steps
    ok my $cnv_pipeline = VRPipe::Pipeline->create(name => 'bcftools_cnv_caller'), 'able to get the bcftools_cnv_caller pipeline';
    my @step_names;
    foreach my $stepmember ($cnv_pipeline->step_members) {
        push(@step_names, $stepmember->step->name);
    }
    is_deeply \@step_names, [qw(vcf_merge_different_samples_control_aware bcftools_cnv_caller)], 'the bcftools_cnv_caller pipeline has the correct steps';
    
    my $cnv_setup = VRPipe::PipelineSetup->create(
        name        => 'cnv_calling',
        pipeline    => $cnv_pipeline,
        datasource  => $vrpipe_ds,
        output_root => $output_dir,
        options     => {}
    );
    
    # figure out output files
    my (@merged_vcfs, @summary_files);
    foreach my $element_id (41 .. 45) {
        my @output_subdirs = output_subdirs($element_id, $cnv_setup->id);
        push(@merged_vcfs, file(@output_subdirs, '1_vcf_merge_different_samples_control_aware', 'merged.vcf.gz'));
        if ($element_id == 41) {
            foreach my $sub_dir (qw(fpdr-fpdr_1 fpdr-fpdr_2 fpdr-fpdr_3)) {
                push(@summary_files, file(@output_subdirs, '2_bcftools_cnv_caller', $sub_dir, 'summary.tab'));
            }
        }
    }
    ok handle_pipeline(@merged_vcfs, @summary_files), 'bcftools_cnv_caller pipeline ran ok and produced the expected output files';
    
    my $results = $summary_files[0]->slurp();
    like $results, qr/AR\t1\t145395605\t145725689\t2\t4\n/, 'the results file for fpdj-fpdj_1 had the expected lines';
    my $cnv_vrfile = VRPipe::File->get(path => $summary_files[0]);
    is $cnv_vrfile->meta_value('sample_control'), 'fpdk_qc1hip5529689', 'sample_control metadata exists on the file';
}

# we'll also take the opportunity to test the loh caller pipeline, since that
# also uses files from the genome studio import
SKIP: {
    my $num_tests = 6;
    skip "hipsci loh calling tests disabled without hipsci_loh_caller.pl in your path", $num_tests unless can_execute('hipsci_loh_caller.pl');
    
    $output_dir = get_output_dir('hipsci_loh_caller');
    
    # check pipeline has correct steps
    ok my $loh_pipeline = VRPipe::Pipeline->create(name => 'hipsci_loh_caller'), 'able to get the hipsci_loh_caller pipeline';
    my @step_names;
    foreach my $stepmember ($loh_pipeline->step_members) {
        push(@step_names, $stepmember->step->name);
    }
    is_deeply \@step_names, [qw(vcf_merge_different_samples_control_aware hipsci_loh_caller)], 'the hipsci_loh_caller pipeline has the correct steps';
    
    my $loh_setup = VRPipe::PipelineSetup->create(
        name        => 'loh_calling',
        pipeline    => $loh_pipeline,
        datasource  => $vrpipe_ds,
        output_root => $output_dir,
        options     => {}
    );
    
    # figure out output files
    my (@merged_vcfs, @loh_files_with_results, @loh_files_no_results);
    foreach my $element_id (41 .. 45) {
        my @output_subdirs = output_subdirs($element_id, $loh_setup->id);
        push(@merged_vcfs, file(@output_subdirs, '1_vcf_merge_different_samples_control_aware', 'merged.vcf.gz'));
        
        my $result_file = file(@output_subdirs, '2_hipsci_loh_caller', 'merged.txt');
        if ($element_id == 41 || $element_id == 45) {
            push(@loh_files_with_results, $result_file);
        }
        else {
            push(@loh_files_no_results, $result_file);
        }
    }
    ok handle_pipeline(@merged_vcfs, @loh_files_with_results), 'hipsci_loh_caller pipeline ran ok and produced the expected output files';
    
    my $created_empty = 0;
    foreach my $file (@loh_files_no_results) {
        $created_empty++ if (-e $file && !-s $file);
    }
    is $created_empty, 3, 'it also produced empty result files for the good cohorts';
    
    my $results  = $loh_files_with_results[0]->slurp();
    my $expected = file(qw(t data fpdk_loh_results.txt))->slurp();
    is $results, $expected, 'the results file for fpdk had the expected lines';
    my $loh_vrfile = VRPipe::File->get(path => $loh_files_with_results[0]);
    is $loh_vrfile->meta_value('sample_control'), 'fpdk_qc1hip5529689', 'sample_control metadata exists on the file';
}

finish;
exit;

sub vcf_lines {
    my $vcf_file = shift;
    my @lines;
    open(my $fh, "zcat $vcf_file |");
    while (<$fh>) {
        next if /^##fileDate/;
        push(@lines, $_);
    }
    close($fh);
    return @lines;
}
