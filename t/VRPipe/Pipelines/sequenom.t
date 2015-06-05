#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 45;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_AUTHOR_TESTS)],
    );
    use TestPipelines;
    use_ok('VRPipe::Schema');
}

ok my $pipeline = VRPipe::Pipeline->create(name => 'sequenom_import_from_irods_and_covert_to_vcf'), 'able to get the sequenom_import_from_irods_and_covert_to_vcf pipeline';
my @step_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@step_names, $stepmember->step->name);
}
is_deeply \@step_names, [qw(irods_get_files_by_basename sequenom_csv_to_vcf)], 'the sequenom_import_from_irods_and_covert_to_vcf pipeline has the correct steps';

my $output_root = get_output_dir('sequenom_import');
my $plex_storage_dir = dir($output_root, 'plex_storage_dir');
mkdir($plex_storage_dir);

my $file_query_file = file(qw(t data datasource.sequenom_fluidigm.filequeries))->absolute->stringify;

ok my $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'seq',
    options => {
        file_query     => $file_query_file,
        local_root_dir => $output_root
    }
  ),
  'could create the irods datasource';

my $ps = VRPipe::PipelineSetup->create(
    name        => 'sequenom & fluidigm import',
    datasource  => $ds,
    output_root => $output_root,
    pipeline    => $pipeline,
    options     => { plex_storage_dir => $plex_storage_dir }
);

my (@vcf_files, @index_files);
my $element_id = 1;
foreach my $basname (qw(QC288261____20130701_G01 QC288261____20130701_C01 QC288261____20130701_A01 QC288261____20130701_E01 S022_1662031438 S023_1662031438 S046_1662031438 S047_1662031438 S070_1662031438 S071_1662031438 S093_1662031438 S094_1662031438 S117_1662031438 S118_1662031438 S141_1662031438 S142_1662031438 S165_1662031438 S166_1662031438 S189_1662031438 S190_1662031438)) {
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

my (@merged_vcf_files, @gtypex_files, @expected_metadata, $gtypex_file);
foreach my $element_id (21 .. 27) {
    my @output_subdirs = output_subdirs($element_id, 2);
    push(@merged_vcf_files, file(@output_subdirs, '1_vcf_merge_different_samples', 'merged.vcf.gz'));
    my $this_gtypex_file = file(@output_subdirs, '2_vcf_genotype_comparison', 'merged.vcf.gz.gtypex');
    push(@gtypex_files, $this_gtypex_file);
    
    my $individual = VRPipe::DataElement->get(id => $element_id)->metadata->{group};
    my %expected = (sample_cohort => $individual);
    if ($individual eq '20f8a331-69ac-4510-94ab-e3a69c50e46f') {
        $expected{genotype_maximum_deviation} = "0.000000:HPSI0813i-ffdb_4_QC1Hip-2";
        $expected{calculated_gender}          = 'M';
        $expected{sample}                     = [qw(QC1Hip-1 QC1Hip-2)];
        $expected{public_name}                = [qw(HPSI0813i-ffdb_3 HPSI0813i-ffdb_4)];
    }
    elsif ($individual eq '3d52354f-8d84-457d-a668-099a758f0e7b') {
        $expected{genotype_maximum_deviation} = '0.000000:HPSI0913i-lofv_33_QC1Hip-4';
        $expected{calculated_gender}          = 'F';
        $expected{sample}                     = 'QC1Hip-4';
        $expected{public_name}                = 'HPSI0913i-lofv_33';
    }
    elsif ($individual eq 'a659498c-81a5-4b04-80ce-36da30632ccf') {
        $expected{genotype_maximum_deviation} = "0.000000:HPSI0813i-ffdc_5_QC1Hip-3";
        $expected{calculated_gender}          = 'M';
        $expected{sample}                     = 'QC1Hip-3';
        $expected{public_name}                = 'HPSI0813i-ffdc_5';
    }
    elsif ($individual eq '693c8457-7595-40bd-a533-22934b02c4b7') {
        $expected{genotype_maximum_deviation} = "5500.000000:HPSI0813pf-piun_QC1Hip-1992";
        $expected{calculated_gender}          = 'F';
        $expected{sample}                     = [qw(QC1Hip-1993 QC1Hip-1995 QC1Hip-1992 QC1Hip-1994)];
        $expected{public_name}                = [qw(HPSI0813i-piun_1 HPSI0813i-piun_3 HPSI0813pf-piun HPSI0813i-piun_2)];
    }
    elsif ($individual eq '9284dc86-5454-473e-b557-dcc783590fa7') {
        $expected{genotype_maximum_deviation} = "0.000000:HPSI0513i-cuau_1_QC1Hip-2001";
        $expected{calculated_gender}          = 'F';
        $expected{sample}                     = [qw(QC1Hip-2000 QC1Hip-2002 QC1Hip-2003 QC1Hip-2001)];
        $expected{public_name}                = [qw(HPSI0513pf-cuau HPSI0513i-cuau_2 HPSI0513i-cuau_3 HPSI0513i-cuau_1)];
        $gtypex_file                          = $this_gtypex_file;
    }
    elsif ($individual eq 'e423ea7c-64b7-43e3-8909-bf290c3846c0') {
        $expected{genotype_maximum_deviation} = "0.000000:HPSI0713i-kaks_2_QC1Hip-1990";
        $expected{calculated_gender}          = 'M';
        $expected{sample}                     = [qw(QC1Hip-1989 QC1Hip-1991 QC1Hip-1988 QC1Hip-1990)];
        $expected{public_name}                = [qw(HPSI0713i-kaks_1 HPSI0713i-kaks_3 HPSI0713pf-kaks HPSI0713i-kaks_2)];
    }
    elsif ($individual eq 'fde0774b-cecd-45d0-9d72-49bd257dd232') {
        $expected{genotype_maximum_deviation} = "0.000000:HPSI0513i-coio_2_QC1Hip-1998";
        $expected{calculated_gender}          = 'F';
        $expected{sample}                     = [qw(QC1Hip-1997 QC1Hip-1999 QC1Hip-1996 QC1Hip-1998)];
        $expected{public_name}                = [qw(HPSI0513i-coio_1 HPSI0513i-coio_3 HPSI0513pf-coio HPSI0513i-coio_2)];
    }
    push(@expected_metadata, \%expected);
}

ok handle_pipeline(@merged_vcf_files, @gtypex_files), 'vcf_merge_and_compare_genotypes pipeline ran ok';

foreach my $vcf_path (@merged_vcf_files) {
    my $meta = VRPipe::File->get(path => $vcf_path)->metadata;
    my $expected = shift @expected_metadata;
    foreach my $key (qw(sample_cohort genotype_maximum_deviation calculated_gender sample public_name)) {
        my $got = ref $meta->{$key}     ? [sort(@{ $meta->{$key} })]    : $meta->{$key};
        my $exp = ref $expected->{$key} ? [sort @{ $expected->{$key} }] : $expected->{$key};
        if (ref $exp) {
            is_deeply $got, $exp, "$key metadata was correct for one of the merged VCF files";
        }
        else {
            is $got, $exp, "$key metadata was correct for one of the merged VCF files";
        }
    }
}

my $schema = VRPipe::Schema->create("VRTrack");
ok my $gtypex_graph_node = $schema->get_file($gtypex_file), 'one of the gtypex files was in the graph db';
my @disc = $gtypex_graph_node->related(outgoing => { type => 'parsed', schema => 'VRTrack', label => 'Discordance' });
is scalar(@disc), 4, 'it had 4 associated Discordance nodes';
my $disc_to_sample = 0;
foreach my $disc (@disc) {
    my ($sample) = $disc->related(incoming => { type => 'discordance', schema => 'VRTrack', label => 'Sample' });
    if ($sample) {
        my $expected = $disc->md5_sample;
        $expected =~ s/^[^_]+_//;
        $disc_to_sample++ if $sample->name eq $expected;
        
        if ($expected eq 'QC1Hip-2000') {
            my $cns_json = $disc->cns;
            my $cns      = $schema->graph->json_decode($cns_json);
            $cns->{type} = $disc->type;
            is_deeply $cns,
              {
                type                                                         => 'fluidigm',
                'HPSI0513i-cuau_2_QC1Hip-2002.HPSI0513pf-cuau_QC1Hip-2000.1' => ['0', '22', '1.00', 'QC1Hip-2002'],
                'HPSI0513i-cuau_3_QC1Hip-2003.HPSI0513pf-cuau_QC1Hip-2000.1' => ['0', '22', '1.00', 'QC1Hip-2003'],
                'HPSI0513i-cuau_1_QC1Hip-2001.HPSI0513pf-cuau_QC1Hip-2000.1' => ['0', '23', '1.00', 'QC1Hip-2001']
              },
              'a discordance node had the expected results stored inside';
        }
    }
}
is $disc_to_sample, 4, 'the 4 discordance nodes were each attached to the correct sample nodes';

finish;
exit;
