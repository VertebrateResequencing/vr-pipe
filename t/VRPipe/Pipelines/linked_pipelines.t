#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 8;
    use VRPipeTest (required_exe => [qw(cat)]);
    use TestPipelines;
}

my $output_dir   = get_output_dir('base_test_pipeline');
my $output_dir_2 = get_output_dir('link_test_pipeline');
my $output_dir_3 = get_output_dir('link_merge_test_pipeline');

ok my $pipeline = VRPipe::Pipeline->get(name => 'test_pipeline'), 'able to get the test_pipeline pipeline';
ok my $link_pipeline = VRPipe::Pipeline->get(name => 'linked_pipeline'), 'able to get the linked_pipeline pipeline';

my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_base_step_names = qw(test_step_one test_step_two test_step_three test_step_four);
is_deeply \@s_names, \@expected_base_step_names, 'the base test pipeline has the correct steps';

my @l_names;
foreach my $stepmember ($link_pipeline->steps) {
    push(@l_names, $stepmember->step->name);
}
my @expected_link_step_names = qw(text_merge test_step_one);
is_deeply \@l_names, \@expected_link_step_names, 'the linked test pipeline has the correct steps';

my $vrobj = VRPipe::Manager->get;
my $tmp_dir = $vrobj->tempdir;

my $orig_ds = file(qw(t data datasource.fofn2));
my $ds = file($tmp_dir, 'datasource.fofn');
copy($orig_ds, $ds);

# setup base pipeline from a file datasource
my $test_pipelinesetup = VRPipe::PipelineSetup->get(name => 'fofn source',
                                                    datasource => VRPipe::DataSource->get(type => 'fofn',
                                                                                          method => 'all',
                                                                                          source => $ds),
                                                    output_root => $output_dir,
                                                    pipeline => $pipeline,
                                                    options => {all_option => 'foo',
                                                                one_option => 50,
                                                                four_option => 'bar',
                                                                cleanup => 0});

# test a vrpipe datasource
my $test_pipelinesetup_2 = VRPipe::PipelineSetup->get(name => 'vrpipe datasource test',
                                                      datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                                            method => 'all',
                                                                                            source => 'fofn source[3]',
                                                                                            options => { }),
                                                      output_root => $output_dir_2,
                                                      pipeline => $link_pipeline,
                                                      options => {all_option => 'bar',
                                                                  one_option => 100,
                                                                  cleanup => 0});

# test a vrpipe datasource, grouping by metadata keys
my $test_pipelinesetup_3 = VRPipe::PipelineSetup->get(name => 'vrpipe datasource grouped test',
                                                      datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                                            method => 'group_by_metadata',
                                                                                            source => '1[2]',
                                                                                            options => { metadata_keys => 'one_meta',
                                                                                                         filter => 'one_meta#50' }),
                                                      output_root => $output_dir_3,
                                                      pipeline => $link_pipeline,
                                                      options => {all_option => 'foo',
                                                                  one_option => 75,
                                                                  cleanup => 0});


# check the three pipelines ran as expected
my (@base_files, @link_files, @link_merge_files);
my $element_id = 0;
foreach my $in ('file.txt', 'file2.txt', 'file3.txt') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    my $step_index = 0;
    foreach my $suffix ('step_one', 'step_one.step_two', 'step_one.step_two.step_three', 'step_one.step_two.step_three.step_four') {
        push(@base_files, file(@output_subdirs, ($step_index+1).'_'.$expected_base_step_names[$step_index], "$in.$suffix"));
        $step_index++;
    }
}

handle_pipeline(@base_files);

# linked pipelines only create their dataelements after the previous pipeline
# has finished, so only now can we figure out remaining expected output file
# paths
$element_id++;
my @output_subdirs = output_subdirs($element_id, 3);
my $step_index = 0;
foreach my $suffix ('txt', 'txt.step_one') {
    push(@link_merge_files, file(@output_subdirs, ($step_index+1).'_'.$expected_link_step_names[$step_index], "merged.$suffix"));
    $step_index++;
}
foreach my $in ('file.txt', 'file2.txt', 'file3.txt') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 2);
    my $step_index = 0;
    foreach my $suffix ('txt', 'txt.step_one') {
        push(@link_files, file(@output_subdirs, ($step_index+1).'_'.$expected_link_step_names[$step_index], "merged.$suffix"));
        $step_index++;
    }
}
ok handle_pipeline(@base_files, @link_files, @link_merge_files), 'pipelines ran and created all expected output files';

# change fofn datasource
my $fh = $ds->openr;
my $tmp_file = file($tmp_dir, 'tmp');
my $tfh = $tmp_file->openw;
<$fh>;
while(<$fh>) {
    print $tfh $_;
}
move($tmp_file, $ds);

# check the three pipelines reset properly and produced new files
ok handle_pipeline(@base_files, @link_files, @link_merge_files), 'pipelines with changed datasource ran ok - original and withdrawn files all exist';
my @new_files;
$element_id++;
@output_subdirs = output_subdirs($element_id, 3);
$step_index = 0;
foreach my $suffix ('txt', 'txt.step_one') {
    push(@new_files, file(@output_subdirs, ($step_index+1).'_'.$expected_link_step_names[$step_index], "merged.$suffix"));
    $step_index++;
}
ok handle_pipeline(@new_files), 'new merge files all exist after changing datasource';

my $element1 = VRPipe::DataElement->get(datasource => 2, result => { paths => [ $base_files[2]->stringify ] });
my $element2 = VRPipe::DataElement->get(datasource => 3, result => { paths => [ map { $_->stringify } @base_files[1,5,9] ], group => '50' });
is_deeply [$element1->withdrawn, $element2->withdrawn], [1,1], 'elements correctly withdrawn after changing datasource';

# # change metadata and detect changes
# VRPipe::File->get(path => $base_files[5]->stringify)->add_metadata({one_meta => '20'}, replace_data => 1);
# ok handle_pipeline(@base_files, @link_files, @link_merge_files, @new_files), 'pipelines ran ok after file metadata changed - original and withdrawn files all exist';
# 
# my @newer_files;
# $element_id++;
# @output_subdirs = output_subdirs($element_id, 3);
# $step_index = 0;
# foreach my $suffix ('txt', 'txt.step_one') {
#     push(@newer_files, file(@output_subdirs, ($step_index+1).'_'.$expected_link_step_names[$step_index], "merged.$suffix"));
#     $step_index++;
# }
# ok handle_pipeline(@newer_files), 'new merge files all exist after metadata change';
# 
# my $element3 = VRPipe::DataElement->get(datasource => 3, result => { paths => [ map { $_->stringify } @base_files[5,9] ], group => '50' });
# is $element3->withdrawn, 1, 'elements correctly withdrawn after metadata change';
# 
# # if filter option one_meta#50 had not been set, a further data element 
# # would have been created. check that we have the expected number of elements
# my $elements = $test_pipelinesetup_3->datasource->elements;
# is scalar(@$elements), 1, 'group_by_metadata/filter option worked correctly - element not created';

finish;