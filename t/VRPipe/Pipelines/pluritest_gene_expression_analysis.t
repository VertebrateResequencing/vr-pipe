#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 10;
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
my @sql = VRPipe::File->create(path => file(qw(t data vrtrack_hipsci_qc1_expression.sql))->absolute)->slurp;
foreach my $sql (@sql) {
    print $mysqlfh $sql;
}
close($mysqlfh);

# setup pipeline
my $output_dir = get_output_dir('import_expression_idat_files');
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
is $results, 4, 'got correct number of idat files from the vrtrack db';

#check pipeline has correct steps
ok my $pipeline = VRPipe::Pipeline->create(name => 'bam_import_from_irods'), 'able to get the idat files via the bam_import_from_irods pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(irods_get_files_by_basename)], 'the pipeline has the correct steps';

# create pipeline setup
VRPipe::PipelineSetup->create(
    name        => 'idat import and qc',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        vrtrack_db     => $ENV{VRPIPE_VRTRACK_TESTDB},
        irods_get_zone => 'archive',
    }
);

#get arrays of output files
my @irods_files;
my @lanes = qw(9252616016_K_Grn 9252616016_I_Grn 9252616016_G_Grn 9252616016_B_Grn);
foreach my $lane (@lanes) {
    push(@irods_files, file($irods_dir, $lane . '.idat'));
}

#run pipeline and check outputs
ok handle_pipeline(@irods_files), 'idat import from irods ran ok and correct idat files found';
my $meta = VRPipe::File->get(path => $irods_files[0])->metadata;
is_deeply $meta,
  {
    'analysis_uuid' => 'cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',
    'bases'         => '0',
    'withdrawn'     => '0',
    'population'    => 'Population',
    'library_tag'   => '9252616016_K',
    'paired'        => '0',
    'reads'         => '0',
    'project'       => 'G0325 [gex] Wellcome Trust Strategic Award application - HIPS',
    'library'       => 'qc1hip5529784_9252616016_K',
    'lane_id'       => '97',
    'individual'    => '27af9a9b-01b2-4cb6-acef-ea52d83e3d26',
    'platform'      => 'SLX',
    'center_name'   => 'SC',
    'sample'        => 'FS11.C',
    'expected_md5'  => '45fe61ebdb49e343558e36aff8dd20b9',
    'study'         => '2625',
    'control'       => 'Stem cell',
    'lane'          => '9252616016_K_Grn',
    'species'       => 'Homo sapiens',
    'insert_size'   => '0',
    'storage_path'  => '/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb'
  },
  'metadata correct for one of the idat files';

#Run merge annotation files pipeline using the output idat files from the import:
$output_dir = get_output_dir('pluritest_gene_expression_analysis');

#check pipeline has correct steps
ok my $merge_pipeline = VRPipe::Pipeline->create(name => 'pluritest_gene_expression_analysis'), 'able to get the pluritest_gene_expression_analysis pipeline';
my @sb_names;
foreach my $stepmember ($merge_pipeline->step_members) {
    push(@sb_names, $stepmember->step->name);
}

is_deeply \@sb_names, [qw(pluritest_annotation_profile_files pluritest_reformat_genome_studio_expression_files pluritest_plot_gene_expression pluritest_vrtrack_update_images)], 'the pluritest_gene_expression_analysis pipeline has the correct steps';

my $r_bin_path       = '/software/vertres/installs/R-2.15/R';
my $pluritest_script = '/software/vertres/scripts/pluriTest_commandLine_vrpipe.r';
my $pluritest_data   = '/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/pluritest.RData';

my $annot_merge = VRPipe::PipelineSetup->create(
    name       => 'expression annot merge',
    pipeline   => $merge_pipeline,
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => 'idat import and qc[1]',
        options => { metadata_keys => 'individual' }
    ),
    options => {
        pluritest_script => $pluritest_script,
        pluritest_data   => $pluritest_data,
        r_bin_path       => $r_bin_path,
        vrtrack_db       => $ENV{VRPIPE_VRTRACK_TESTDB}
    },
    output_root => $output_dir,
);

my @final_files;
foreach my $element (@{ get_elements($annot_merge->datasource) }) {
    my @output_dirs = output_subdirs($element->id, $annot_merge->id);
    foreach my $kind (qw(01 02a 02 03c 03)) {
        push(@final_files, file(@output_dirs, '3_pluritest_plot_gene_expression', 'pluritest_image' . $kind . '.png'));
    }
    push(@final_files, file(@output_dirs, '3_pluritest_plot_gene_expression', 'pluritest.csv'));
}
ok handle_pipeline(@final_files), 'pluritest_gene_expression_analysis pipeline ran ok and produced the expected image files';

finish;
