#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 5;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB VRPIPE_AUTHOR_TESTS WAREHOUSE_DATABASE WAREHOUSE_HOST WAREHOUSE_PORT WAREHOUSE_USER VRPIPE_PLURITEST_RSCRIPT VRPIPE_PLURITEST_RDATA VRPIPE_PLURITEST_RBIN VRPIPE_PLURITEST_RLIBS)],
        required_exe => [qw(iget iquest)]
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

my $output_dir = get_output_dir('import_expression_idat_files');
my $irods_dir = dir($output_dir, 'irods_import')->stringify;

# create a setup to import idat files from irods
my $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'archive',
    options => {
        file_query     => q[study_id = 2625 and type = idat and dcterms:created '<' 2013-06-07],
        local_root_dir => $irods_dir
    }
);

my $pipeline = VRPipe::Pipeline->create(name => 'vrtrack_populate_from_irods_and_download_files');

VRPipe::PipelineSetup->create(
    name        => 'idat populate',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        vrtrack_db         => $ENV{VRPIPE_VRTRACK_TESTDB},
        vrlane_storage_dir => $irods_dir
    }
);

my @analysis_files = (file($irods_dir, '/archive/GAPI/exp/analysis/69/61/88/hipsci_7samples_2013-05-15/hipsci_2013-05-15_Sample_Probe_Profile.txt.gz'), file($irods_dir, '/archive/GAPI/exp/analysis/69/61/88/hipsci_7samples_2013-05-15/hipsci_2013-05-15_annotation.txt'));
ok handle_pipeline(@analysis_files), 'vrtrack_populate_from_vrpipe_metadata pipeline ran ok and got the analysis files';

# create pluritest setup using the output idat files from the import
$output_dir = get_output_dir('pluritest_gene_expression_analysis');

# check pipeline has correct steps
ok my $pluritest_pipeline = VRPipe::Pipeline->create(name => 'pluritest_gene_expression_analysis'), 'able to get the pluritest_gene_expression_analysis pipeline';
my @sb_names;
foreach my $stepmember ($pluritest_pipeline->step_members) {
    push(@sb_names, $stepmember->step->name);
}

is_deeply \@sb_names, [qw(pluritest_annotation_profile_files pluritest_reformat_genome_studio_expression_files pluritest_plot_gene_expression pluritest_vrtrack_update_images)], 'the pluritest_gene_expression_analysis pipeline has the correct steps';

my $r_bin_path       = $ENV{VRPIPE_PLURITEST_RBIN};
my $r_libs           = $ENV{VRPIPE_PLURITEST_RLIBS};
my $pluritest_script = $ENV{VRPIPE_PLURITEST_RSCRIPT};
my $pluritest_data   = $ENV{VRPIPE_PLURITEST_RDATA};

my $pluri_setup = VRPipe::PipelineSetup->create(
    name       => 'pluri setup',
    pipeline   => $pluritest_pipeline,
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => 'idat populate[2:input_files]',
        options => { metadata_keys => 'sample_cohort|beadchip' }
    ),
    options => {
        pluritest_script => $pluritest_script,
        pluritest_data   => $pluritest_data,
        r_bin_path       => $r_bin_path,
        r_libs           => $r_libs,
        vrtrack_db       => $ENV{VRPIPE_VRTRACK_TESTDB}
    },
    output_root => $output_dir
);

my @final_files;
foreach my $element (@{ get_elements($pluri_setup->datasource) }) {
    my @output_dirs = output_subdirs($element->id, $pluri_setup->id);
    foreach my $kind (qw(01 02a 02 03c 03)) {
        push(@final_files, file(@output_dirs, '3_pluritest_plot_gene_expression', 'pluritest_image' . $kind . '.png'));
    }
    push(@final_files, file(@output_dirs, '3_pluritest_plot_gene_expression', 'pluritest.csv'));
}
ok handle_pipeline(@final_files), 'pluritest_gene_expression_analysis pipeline ran ok and produced the expected image files';

finish;
exit;
