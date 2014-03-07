#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;
use Data::Dumper;

BEGIN {
    use Test::Most tests => 28;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB)],
        required_exe => [qw(iget iquest fcr-to-vcf sort bgzip)]
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

# (for testing purposes later on in this script, the vrtrack db here has 8
#  samples, 4 per individual, but the library name for lane 9300870057_R02C02
#  (sample fpdr_3, individual 6d3d2acf-29a5-41a2-8992-1414706a527d) was changed
#  from 283163_H01_qc1hip5533833 to 283163_G01_qc1hip5529689 (for an unrelated
#  individual))

# setup pipeline
my $output_dir = get_output_dir('genome_studio_import_from_irods_and_convert_to_vcf');
my $irods_dir = dir($output_dir, 'irods_import')->stringify;

# setup vrtrack datasource
ok my $vrtrack_ds = VRPipe::DataSource->create(
    type    => 'vrtrack',
    method  => 'analysis_genome_studio',
    source  => $ENV{VRPIPE_VRTRACK_TESTDB},
    options => { local_root_dir => $irods_dir }
  ),
  'could create a vrtrack datasource';

# check correct number of gtc file retrieved
my $results = 0;
foreach my $element (@{ get_elements($vrtrack_ds) }) {
    $results++;
}
is $results, 8, 'got correct number of gtc files from the vrtrack db';

# check pipeline has correct steps
ok my $import_pipeline = VRPipe::Pipeline->create(name => 'genome_studio_import_from_irods_and_convert_to_vcf'), 'able to get the genome_studio_import_from_irods_and_convert_to_vcf pipeline';
my @s_names;
foreach my $stepmember ($import_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(irods_get_files_by_basename split_genome_studio_genotype_files illumina_coreexome_manifest_to_map genome_studio_fcr_to_vcf)], 'the pipeline has the correct steps';

my $external_gzip_file        = file(qw(t data hipsci_genotyping.fcr.txt.gz))->absolute->stringify;
my $external_reheader_penncnv = file(qw(t data reheader_penncnv.txt))->absolute->stringify;
my $manifest_file_soure       = file(qw(t data hipsci_genotyping.snp.manifest))->absolute->stringify;
my $manifest_file             = file($output_dir, 'hipsci_genotyping.manifest')->stringify;
copy($manifest_file_soure, $manifest_file);

# create pipeline setup
VRPipe::PipelineSetup->create(
    name        => 'gtc import',
    datasource  => $vrtrack_ds,
    output_root => $output_dir,
    pipeline    => $import_pipeline,
    options     => {
        vrtrack_db         => $ENV{VRPIPE_VRTRACK_TESTDB},
        irods_get_zone     => 'archive',
        external_gzip_file => $external_gzip_file,
        reheader_penncnv   => $external_reheader_penncnv,
        coreexome_manifest => $manifest_file,
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
my @samples    = qw(fpdj fpdj_2 fpdj_1 fpdr fpdr_2 fpdr_1 fpdj_3 fpdr_3);
foreach my $sample (@samples) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@genotype_files, file(@output_subdirs, '2_split_genome_studio_genotype_files', $sample . '.genotyping.fcr.txt'));
    push(@vcf_files,      file(@output_subdirs, '4_genome_studio_fcr_to_vcf',           "$sample/$sample.vcf.gz"));
}

# run pipeline and check outputs
ok handle_pipeline(@genotype_files, @vcf_files), 'bam import from irods and split genome studio genotype file pipeline ran ok';

# check genotype file metadata
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
is_deeply [vcf_lines($vcf_files[0])], [vcf_lines(file(qw(t data hipsci_genotyping.fpdj.snp.vcf.gz)))], 'VCF file produced was as expected';
$meta = VRPipe::File->get(path => $vcf_files[0])->metadata;
is $meta->{individual}, 'ca04b23b-c5b0-4389-95a3-5c7c8e6d51f2', 'the VCF file has individual metadata';

# we'll take this opportunity to test the vcf_merge_and_compare_genotypes
# pipeline as well
$output_dir = get_output_dir('vcf_merge_and_compare_genotypes');
my $vrpipe_ds = VRPipe::DataSource->create(
    type    => 'vrpipe',
    method  => 'group_by_metadata',
    source  => '1[4]',
    options => { metadata_keys => 'individual' }
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
foreach my $element_id (9, 10) {
    my @output_subdirs = output_subdirs($element_id, 2);
    push(@genotype_files,   file(@output_subdirs, '2_vcf_merge_different_samples', 'merged.vcf.gz'));
    push(@merged_vcf_files, file(@output_subdirs, '4_vcf_genotype_comparison',     'merged.vcf.gz.gtypex'));
    
    my $individual = VRPipe::DataElement->get(id => $element_id)->metadata->{group};
    my %expected = (individual => $individual);
    if ($individual eq 'ca04b23b-c5b0-4389-95a3-5c7c8e6d51f2') {
        $expected{genotype_maximum_deviation} = ['==', 0, 'fpdj_3'];
    }
    else {
        $expected{genotype_maximum_deviation} = ['>', 10, 'fpdr_3'];
    }
    push(@expected_metadata, \%expected);
}

ok handle_pipeline(@merged_vcf_files, @gtypex_files), 'vcf_merge_and_compare_genotypes pipeline ran ok';

foreach my $vcf_path (@merged_vcf_files) {
    $meta = VRPipe::File->get(path => $vcf_path)->metadata;
    my $expected = shift @expected_metadata;
    
    is $meta->{individual}, $expected->{individual}, "individual metadata was correct for one of the merged VCF files";
    
    my ($cmp, $eval, $esample) = @{ $expected->{genotype_maximum_deviation} };
    my ($aval, $asample) = split(':', $meta->{genotype_maximum_deviation});
    cmp_ok $aval, $cmp, $eval, "genotype_maximum_deviation metadata value was correct for one of the merged VCF files";
    is $asample, $esample, "genotype_maximum_deviation metadata sample was correct for one of the merged VCF files";
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
    
    # get array of output files and check outputs as the pipeline is run
    my (@merged_vcfs, $loh_file_with_results, $loh_file_empty);
    foreach my $element_id (9, 10) {
        my @output_subdirs = output_subdirs($element_id, $loh_setup->id);
        push(@merged_vcfs, file(@output_subdirs, '1_vcf_merge_different_samples_control_aware', 'merged.vcf.gz'));
        if ($element_id == 9) {
            $loh_file_with_results = file(@output_subdirs, '2_hipsci_loh_caller', 'merged.txt');
        }
        else {
            $loh_file_empty = file(@output_subdirs, '2_hipsci_loh_caller', 'merged.txt');
        }
    }
    ok handle_pipeline(@merged_vcfs, $loh_file_with_results), 'hipsci_loh_caller pipeline ran ok and produced the expected output files';
    
    my $empty_created = -e $loh_file_empty && !-s $loh_file_empty;
    ok $empty_created, 'it also produced an empty result file for fpdj';
    
    my $results = $loh_file_with_results->slurp();
    is $results, "1\t9304731\t29355519\tfpdr_3\t1956 (SNPs)\n6\t29496664\t78649193\tfpdr_3\t1019 (SNPs)\n", 'the results file for fpdr had the expected lines';
    my $loh_vrfile = VRPipe::File->get(path => $loh_file_with_results);
    is $loh_vrfile->meta_value('control'), 'fpdr', 'control metadata exists on the file';
}

# we'll also take the opportunity to test the penncnv pipeline, since that also
# uses files from the genome studio import
SKIP: {
    my $num_tests = 4;
    skip "penncnv cnv calling tests disabled without VRPIPE_PENNCNV_SCRIPT_DIR and VRPIPE_PENNCNV_LIB_DIR VRPIPE_PENNCNV_PERL environment variables", $num_tests unless ($ENV{VRPIPE_PENNCNV_SCRIPT_DIR} && $ENV{VRPIPE_PENNCNV_LIB_DIR} && $ENV{VRPIPE_PENNCNV_PERL});
    
    $output_dir = get_output_dir('penncnv_cnv_calling');
    
    # check pipeline has correct steps
    ok my $penn_pipeline = VRPipe::Pipeline->create(name => 'penncnv_cnv_calling'), 'able to get the penncnv_cnv_calling pipeline';
    my @sp_names;
    foreach my $stepmember ($penn_pipeline->step_members) {
        push(@sp_names, $stepmember->step->name);
    }
    is_deeply \@sp_names, [qw(penncnv_detect_cnv penncnv_filter_cnv)], 'the penncnv_cnv_calling pipeline has the correct steps';
    
    # so that these tests will pass, we need to use a new gtc import based on
    # a different bit of the external_gzip_file (which is too big to use all of,
    # and where this new bit fails the previous tests)
    $external_gzip_file = file(qw(t data hipsci_genotyping.fcr_v2.txt.gz))->absolute->stringify;
    my $gtc_import_setup = VRPipe::PipelineSetup->create(
        name        => 'gtc import 2',
        datasource  => $vrtrack_ds,
        output_root => $output_dir,
        pipeline    => $import_pipeline,
        options     => {
            vrtrack_db         => $ENV{VRPIPE_VRTRACK_TESTDB},
            irods_get_zone     => 'archive',
            external_gzip_file => $external_gzip_file,
            reheader_penncnv   => $external_reheader_penncnv,
            coreexome_manifest => $manifest_file,
            cleanup            => 0
        }
    );
    handle_pipeline();
    
    # create pipeline setup
    my $penncnv_script_dir  = $ENV{VRPIPE_PENNCNV_SCRIPT_DIR};
    my $detect_cnv_script   = file($penncnv_script_dir, 'detect_cnv.pl');
    my $filter_cnv_script   = file($penncnv_script_dir, 'filter_cnv.pl');
    my $penncnv_lib_dir     = $ENV{VRPIPE_PENNCNV_LIB_DIR};
    my $penncnv_perl        = $ENV{VRPIPE_PENNCNV_PERL};
    my $detect_cnv_hmm      = file($penncnv_lib_dir, 'custom.hmm');
    my $detect_cnv_pfb      = file($penncnv_lib_dir, 'HumanExome12v1.1.hg19.pfb');
    my $filter_numsnps      = 8;
    my $filter_length       = '120k';
    my $filter_confidence   = 8;
    my $gtc_import_setup_id = $gtc_import_setup->id;
    my $penn_setup          = VRPipe::PipelineSetup->create(
        name       => 'penncnv_calling',
        pipeline   => $penn_pipeline,
        datasource => VRPipe::DataSource->create(
            type   => 'vrpipe',
            method => 'all',
            source => $gtc_import_setup_id . '[2]',
        ),
        output_root => $output_dir,
        options     => {
            detect_cnv_script    => $detect_cnv_script,
            detect_cnv_hmm       => $detect_cnv_hmm,
            detect_cnv_pfb       => $detect_cnv_pfb,
            filter_cnv_script    => $filter_cnv_script,
            filter_numsnps       => $filter_numsnps,
            filter_length        => $filter_length,
            filter_confidence    => $filter_confidence,
            perl_for_penncnv_exe => $penncnv_perl,
        }
    );
    
    # get array of output files and check outputs as the pipeline is run
    my @penncnv_files;
    my $pager = $penn_setup->datasource->elements;
    my @element_ids;
    while (my $es = $pager->next) {
        push(@element_ids, map { $_->id } @$es);
    }
    
    foreach my $sample (@samples) {
        my @output_subdirs = output_subdirs(shift(@element_ids), $penn_setup->id);
        push(@penncnv_files, file(@output_subdirs, '2_penncnv_filter_cnv', $sample . '.genotyping.fcr.txt.rawcnv.filtercnv'));
    }
    ok handle_pipeline(@penncnv_files), 'penncnv_cnv_calling pipeline ran ok and produced the expected output files';
    
    # check cnv file metadata
    my $penn_meta = VRPipe::File->get(path => $penncnv_files[0])->metadata;
    is_deeply $penn_meta, {
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
        
        'sample'       => 'fpdj',
        'expected_md5' => 'a965097662c37a2d74a7ca429f58f726',
        'study'        => '2624',
        'control'      => 'Control',
        'lane'         => '9300870057_R06C02',
        
        'species'     => 'Homo sapiens',
        'insert_size' => '0',
        
        'storage_path' => '/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',
        
        'cnv_analysis_type' => 'penncnv'
      
      },
      'metadata correct for one of the penncnv files';
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
