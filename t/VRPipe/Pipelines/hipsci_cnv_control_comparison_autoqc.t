#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Cwd;
use Data::Dumper;

BEGIN {
    use Test::Most tests => 9;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB)],
    );
    use TestPipelines;
    use_ok('VRTrack::Factory');
}

#create database for autoqc to update:
my %cd = VRTrack::Factory->connection_details('rw');
open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
my @sql = VRPipe::File->create(path => file(qw(t data vrtrack_hipsci_qc1_genotyping_autoqc.sql))->absolute)->slurp;
foreach my $sql (@sql) {
    print $mysqlfh $sql;
}
close($mysqlfh);

#Run quantisnp pipeline using the output genotype files from the import:
my $output_dir = get_output_dir('cnv_control_comparison');

#check pipeline has correct steps
ok my $bed_pipeline = VRPipe::Pipeline->create(name => 'cnv_control_comparison'), 'able to get the cnv_control_comparison pipeline';
my @sb_names;
foreach my $stepmember ($bed_pipeline->step_members) {
    push(@sb_names, $stepmember->step->name);
}
is_deeply \@sb_names, [qw(reformat_cnv_output_to_bed cnv_control_comparison)], 'the cnv_control_comparison pipeline has the correct steps';

ok my $ds = VRPipe::DataSource->create(
    type   => 'fofn',
    method => 'group_all',
    source => file(qw(t data hipsci_control_removal.fofn))->absolute->stringify
  ),
  'could create a fofn datasource';

my $suffix = ".filtercnv";
my $cohort = '2a39941c-12b2-41bf-92f3-70b88b66a3a4';
my $uuid   = '12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab';

my %sample_lanes = (
    'qc1hip5533830' => '9300870166_R06C01',
    'qc1hip5533829' => '9300870166_R05C01',
    'qc1hip5533832' => '9300870166_R02C02',
    'qc1hip5533831' => '9300870166_R01C02'
);

foreach my $sample (qw(qc1hip5533830 qc1hip5533831 qc1hip5533829 qc1hip5533832)) {
    my $file = VRPipe::File->create(path => file('t', 'data', $sample . $suffix)->absolute);
    my $control = $sample eq 'qc1hip5533832' ? 'Control' : 'Stem cell';
    $file->add_metadata({ sample        => $sample });
    $file->add_metadata({ control       => $control });
    $file->add_metadata({ individual    => $cohort });
    $file->add_metadata({ analysis_uuid => $uuid });
    $file->add_metadata({ lane          => $sample_lanes{$sample} });
}

my $quan_bed = VRPipe::PipelineSetup->create(
    name        => 'penncnv_reformat_bed',
    pipeline    => $bed_pipeline,
    datasource  => $ds,
    output_root => $output_dir,
    options     => {
        cnv_analysis_type => 'penncnv',
    }
);

my @pennbed_files;
my $element_id = 1;
foreach my $sample (qw(qc1hip5533830 qc1hip5533829 qc1hip5533831)) {
    my @output_subdirs = output_subdirs($element_id, 1);
    push(@pennbed_files, file(@output_subdirs, '2_cnv_control_comparison', $cohort . '_' . $sample . '_penncnv.bed.intersection'));
    push(@pennbed_files, file(@output_subdirs, '2_cnv_control_comparison', $cohort . '_' . $sample . '_penncnv.bed.diff'));
}
ok handle_pipeline(@pennbed_files), 'cnv_control_comparison pipeline ran ok and produced the expected output files';

#check cnv file metadata
my $meta = VRPipe::File->get(path => $pennbed_files[1])->metadata;
is_deeply $meta,
  {
    'analysis_uuid'     => '12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',
    'cnv_total'         => '25',
    'cnv_minus_control' => '20',
    'sample'            => 'qc1hip5533830',
    'control'           => 'Stem cell',
    'individual'        => '2a39941c-12b2-41bf-92f3-70b88b66a3a4',
    'lane'              => '9300870166_R06C01'
  },
  'intersection and diff metadata successfully added to one of the penncnv bed files';

#Run autoqc pipeline using the output bed files from the cnv_control_removal:
$output_dir = get_output_dir('vrtrack_update_mapstats_hipsci');

#check pipeline has correct steps
ok my $auto_pipeline = VRPipe::Pipeline->create(name => 'vrtrack_update_mapstats_hipsci'), 'able to get the vrtrack_update_mapstats_hipsci pipeline';
my @sq_names;
foreach my $stepmember ($auto_pipeline->step_members) {
    push(@sq_names, $stepmember->step->name);
}
is_deeply \@sq_names, [qw(vrtrack_update_mapstats_hipsci)], 'the vrtrack_update_mapstats_hipsci pipeline has the correct steps';

my $penn_auto = VRPipe::PipelineSetup->create(
    name       => 'penncnv_auto_qc',
    pipeline   => $auto_pipeline,
    datasource => VRPipe::DataSource->create(
        type   => 'vrpipe',
        method => 'all',
        source => 'penncnv_reformat_bed[2]'
    ),
    output_root => $output_dir,
    options     => { vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB} }
);

#Get array of output files and check outputs as the pipeline is run
ok handle_pipeline(), 'vrtrack_update_mapstats_hipsci pipeline ran ok';

done_testing;
exit;
