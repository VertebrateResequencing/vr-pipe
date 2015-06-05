#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 13;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_AUTHOR_TESTS WAREHOUSE_DATABASE WAREHOUSE_HOST WAREHOUSE_PORT WAREHOUSE_USER VRPIPE_PLURITEST_RSCRIPT VRPIPE_PLURITEST_RDATA VRPIPE_PLURITEST_RBIN VRPIPE_PLURITEST_RLIBS)],
        required_exe => [qw(iget iquest)]
    );
    use TestPipelines;
    use_ok('VRPipe::Schema');
}

my $output_dir = get_output_dir('import_expression_idat_files');
my $irods_dir = dir($output_dir, 'irods_import')->stringify;

# create a setup to import idat files from irods
my $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'archive',
    options => {
        file_query     => q[study_id = 2625 and type = idat and sample_cohort = fee2bf14-59cd-4c57-b70e-12941437b0d5],
        local_root_dir => $irods_dir
    }
);

my $pipeline = VRPipe::Pipeline->create(name => 'irods_analysis_files_download');

VRPipe::PipelineSetup->create(
    name        => 'idat import',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        irods_download_input_files => 1,
        vrlane_storage_dir         => $irods_dir
    }
);

my @analysis_files = (file($irods_dir, '/archive/GAPI/exp/analysis/be/9a/98/hipsci_21samples_2013-08-30/hipsci_21samples_2013-08-30_Sample_Probe_Profile.txt'), file($irods_dir, '/archive/GAPI/exp/analysis/be/9a/98/hipsci_21samples_2013-08-30/hipsci_21samples_2013-08-30_annotation.txt'));
ok handle_pipeline(@analysis_files), 'irods_analysis_files_download pipeline ran ok and got the analysis files';

# create pluritest setup using the output idat files from the import
$output_dir = get_output_dir('pluritest_gene_expression_analysis');

# check pipeline has correct steps
ok my $pluritest_pipeline = VRPipe::Pipeline->create(name => 'pluritest_gene_expression_analysis'), 'able to get the pluritest_gene_expression_analysis pipeline';
my @sb_names;
foreach my $stepmember ($pluritest_pipeline->step_members) {
    push(@sb_names, $stepmember->step->name);
}

is_deeply \@sb_names, [qw(pluritest_annotation_profile_files pluritest_reformat_genome_studio_expression_files pluritest_plot_gene_expression)], 'the pluritest_gene_expression_analysis pipeline has the correct steps';

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
        source  => 'idat import[1:input_files]',
        options => { metadata_keys => 'sample_cohort' }
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

my (@plot_files, @csv_files);
foreach my $element (@{ get_elements($pluri_setup->datasource) }) {
    my @output_dirs = output_subdirs($element->id, $pluri_setup->id);
    foreach my $kind (qw(01 02a 02 03c 03)) {
        push(@plot_files, file(@output_dirs, '3_pluritest_plot_gene_expression', 'pluritest_image' . $kind . '.png'));
    }
    push(@csv_files, file(@output_dirs, '3_pluritest_plot_gene_expression', 'pluritest.csv'));
}
ok handle_pipeline(@plot_files, @csv_files), 'pluritest_gene_expression_analysis pipeline ran ok and produced the expected image files';

my $vrtrack = VRPipe::Schema->create("VRTrack");
ok my $plot_node = $vrtrack->get_file($plot_files[0]), 'one of the plots was in the graph db';
is $plot_node->property('type'), 'intensity', 'the plot file had the correct type';
my ($donor_node) = $plot_node->related(incoming => { namespace => 'VRTrack', label => 'Donor', type => 'pluritest_plot' });
is $donor_node->id, 'fee2bf14-59cd-4c57-b70e-12941437b0d5', 'the plot was related to the correct Donor node';
my (@plots) = $donor_node->related(outgoing => { type => 'pluritest_plot' });
is scalar(@plots), 5, 'all 5 plots were attached to the donor node';

ok my $csv_node = $vrtrack->get_file($csv_files[0]), 'the csv was in the graph db';
my @plur_nodes = $csv_node->related(outgoing => { namespace => 'VRTrack', label => 'Pluritest', type => 'parsed' });
is scalar(@plur_nodes), 7, 'the csv was attached to 7 Pluritest nodes';
my ($plur_node) = grep { $_->md5_sample =~ /QC1Hip-86/ } @plur_nodes;
my ($sample_node) = $plur_node->related(incoming => { namespace => 'VRTrack', label => 'Sample', type => 'pluritest' });
is $sample_node->name, 'QC1Hip-86', 'a Pluritest node was attached to the correct sample';

my $data = $vrtrack->graph->json_decode($plur_node->data);
is_deeply $data, { 'pluri-raw' => '20.435', 'pluri logit-p' => 1, novelty => '1.647', 'novelty logit-p' => '0.103', RMSD => '0.757' }, 'a Pluritest node had the correct data';

finish;
exit;
