#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Cwd;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('quantisnp_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'quantisnp'), 'able to get the quantisnp pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(quantisnp_reformat_gs_export quantisnp_detect_cnv);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

ok my $ds = VRPipe::DataSource->create(
    type   => 'fofn',
    method => 'all',
    source => file(qw(t data hipsci_cnv_calling.fofn))->absolute
  ),
  'could create a fofn datasource';

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name        => 'my quantisnp_pipeline',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        #'reformat_annotation' => file(qw(t data hipsci_gene_expression_annotation.txt))->absolute->stringify,
        #'reformat_mapping'    => file(qw(t data hipsci_gene_expression_mapping.txt))->absolute->stringify,
        #'pluritest_script'    => $pluritest_script,
        #'pluritest_data'      => file(qw(t data pluritest.RData))->absolute->stringify,
        #'r_bin_path'          => $r_bin_path,
        #'test_bin_path'          => $test_bin_path,
        cleanup => 0
    }
);

my (@output_files, @final_files);
my @output_dirs = output_subdirs(1);
#push(@output_files, file(@output_dirs, '1_penncnv_detect_cnv', "hipsci_cnv_calling.txt.rawcnv"));
#foreach my $kind (qw(01 02a 02 03c 03)) {
#    push(@final_files, file(@output_dirs, '2_t2', 'pluritest_image' . $kind . '.png'));
#}
#push(@final_files, file(@output_dirs, '2_penncnv_filter_cnv', 'hipsci_cnv_calling.txt.filtercnv'));

print "########## QUANTISNP - CALLING CNVS #################\n";

#print "Expected output files: @output_files\n";
#print "Expected final files: @final_files\n";

$#output_files = -1;
$#final_files  = -1;
#push(@output_files, '271298_A01_hipscigt5466706.txt.rawcnv');
push(@output_files, file(@output_dirs, '1_quantisnp_reformat_gs_export', "271298_A01_hipscigt5466706_50000.txt.reformat"));
push(@final_files,  file(@output_dirs, '2_quantisnp_detect_cnv',         "271298_A01_hipscigt5466706_50000.cnv"));

print "Changed output to: @output_files\n";
print "Changed final to: @final_files\n";

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
