#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 21;
    use VRPipeTest (required_exe => [qw(cat)]);
    use TestPipelines;
}

my $output_dir   = get_output_dir('base_test_pipeline');
my $output_dir_2 = get_output_dir('link_test_pipeline');
my $output_dir_3 = get_output_dir('link_merge_test_pipeline');

ok my $pipeline      = VRPipe::Pipeline->create(name => 'test_pipeline'),   'able to get the test_pipeline pipeline';
ok my $link_pipeline = VRPipe::Pipeline->create(name => 'linked_pipeline'), 'able to get the linked_pipeline pipeline';

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

my $vrobj   = VRPipe::Manager->get;
my $tmp_dir = $vrobj->tempdir;

my $orig_ds = file(qw(t data datasource.fofn2));
my $ds = file($tmp_dir, 'datasource.fofn');
copy($orig_ds, $ds);

# setup base pipeline from a file datasource
my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'fofn source',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => $ds
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        all_option  => 'foo',
        one_option  => 50,
        four_option => 'bar',
        cleanup     => 0
    }
);

# test a vrpipe datasource
my $test_pipelinesetup_2 = VRPipe::PipelineSetup->create(
    name       => 'vrpipe datasource test',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'fofn source[3]',
        options => {}
    ),
    output_root => $output_dir_2,
    pipeline    => $link_pipeline,
    options     => {
        all_option => 'bar',
        one_option => 100,
        cleanup    => 0
    }
);

# test a vrpipe datasource, grouping by metadata keys
my $test_pipelinesetup_3 = VRPipe::PipelineSetup->create(
    name       => 'vrpipe datasource grouped test',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => '1[2]',
        options => {
            metadata_keys => 'one_meta',
            filter        => 'one_meta#50'
        }
    ),
    output_root => $output_dir_3,
    pipeline    => $link_pipeline,
    options     => {
        all_option => 'foo',
        one_option => 75,
        cleanup    => 0
    }
);

# check the three pipelines ran as expected
my (@base_files, @link_files, @link_merge_files, @step_four_base_files);
my $element_id = 0;
foreach my $in ('file.txt', 'file2.txt', 'file3.txt') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    my $step_index     = 0;
    foreach my $suffix ('step_one', 'step_one.step_two', 'step_one.step_two.step_three', 'step_one.step_two.step_three.step_four') {
        my $file = file(@output_subdirs, ($step_index + 1) . '_' . $expected_base_step_names[$step_index], "$in.$suffix");
        push(@base_files, $file);
        push(@step_four_base_files, $file) if ($suffix eq 'step_one.step_two.step_three.step_four');
        $step_index++;
    }
}

handle_pipeline(@base_files);
handle_pipeline(); # *** there's some sort of bug that means the following test can fail with only 2 elements sometimes; does this somehow help?

# linked pipelines only create their dataelements after the previous pipeline
# has finished, so only now can we figure out remaining expected output file
# paths
my $element_count = 0;
foreach my $eid (map { $_->id } @{ get_elements($test_pipelinesetup_2->datasource) }) {
    $element_count++;
    my @output_subdirs = output_subdirs($eid, 2);
    my $step_index = 0;
    foreach my $suffix ('txt', 'txt.step_one') {
        push(@link_files, file(@output_subdirs, ($step_index + 1) . '_' . $expected_link_step_names[$step_index], "merged.$suffix"));
        $step_index++;
    }
}
is $element_count, 3, 'the "all" linked pipeline had a dataelement for each parent element';
$element_count = 0;
foreach my $eid (map { $_->id } @{ get_elements($test_pipelinesetup_3->datasource) }) {
    $element_count++;
    my @output_subdirs = output_subdirs($eid, 3);
    my $step_index = 0;
    foreach my $suffix ('txt', 'txt.step_one') {
        push(@link_merge_files, file(@output_subdirs, ($step_index + 1) . '_' . $expected_link_step_names[$step_index], "merged.$suffix"));
        $step_index++;
    }
}
is $element_count, 1, 'the "group_by_metadata" linked pipeline had just 1 element for all the input elements';
ok handle_pipeline(@base_files, @link_files, @link_merge_files), 'pipelines ran and created all expected output files';

# Check that metadata filter is applied after grouping.
# Only one file has the metadata three_meta => StepOption_default_decided_three_option -
# only this one file should be required to pass the filter.
my $ds_test = VRPipe::DataSource->create(
    type    => 'vrpipe',
    method  => 'group_by_metadata',
    source  => '1[4]',
    options => {
        metadata_keys => 'four_meta',
        filter        => 'three_meta#StepOption_default_decided_three_option'
    }
);
my @results = ();
foreach my $element (@{ get_elements($ds_test) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => \@step_four_base_files, group => 'bar' }], 'metadata filtering for "group_by_metadata" vrpipe datasource method worked as expected';

$ds_test = VRPipe::DataSource->create(
    type    => 'vrpipe',
    method  => 'group_by_metadata',
    source  => '1[4]',
    options => {
        metadata_keys         => 'four_meta',
        filter                => 'three_meta#StepOption_default_decided_three_option',
        filter_after_grouping => 0
    }
);
@results = ();
foreach my $element (@{ get_elements($ds_test) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => [$step_four_base_files[1]], group => 'bar' }], 'metadata filtering for "group_by_metadata" vrpipe datasource method with "filter_after_grouping" option off worked as expected';

$ds_test = VRPipe::DataSource->create(
    type    => 'vrpipe',
    method  => 'all',
    source  => '1[4]',
    options => { filter => 'three_meta#StepOption_default_decided_three_option' }
);
@results = ();
foreach my $element (@{ get_elements($ds_test) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => [$step_four_base_files[1]] }], 'metadata filtering for "all" vrpipe datasource method worked as expected';

$ds_test = VRPipe::DataSource->create(
    type    => 'vrpipe',
    method  => 'all',
    source  => '1[0]',
    options => {}
);
@results = ();
foreach my $element (@{ get_elements($ds_test) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => [file(qw(t data file.txt))->absolute->stringify] }, { paths => [file(qw(t data file2.txt))->absolute->stringify] }, { paths => [file(qw(t data file3.txt))->absolute->stringify] }], 'getting vrpipe datasource from step 0 worked';

$ds_test = VRPipe::DataSource->create(
    type    => 'vrpipe',
    method  => 'group_all',
    source  => '1[4]',
    options => {}
);
@results = ();
foreach my $element (@{ get_elements($ds_test) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => \@step_four_base_files }], 'vrpipe datasource group_all method worked as expected';

$ds_test = VRPipe::DataSource->create(
    type    => 'vrpipe',
    method  => 'group_all',
    source  => '1[4]',
    options => { filter => 'three_meta#StepOption_default_decided_three_option' }
);
@results = ();
foreach my $element (@{ get_elements($ds_test) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => [$step_four_base_files[1]] }], 'vrpipe datasource group_all with filter option worked as expected';

# change fofn datasource
my $fh       = $ds->openr;
my $tmp_file = file($tmp_dir, 'tmp');
my $tfh      = $tmp_file->openw;
<$fh>;
while (<$fh>) {
    print $tfh $_;
}
move($tmp_file, $ds);

# check the three pipelines reset properly and produced new files
ok handle_pipeline(@base_files, @link_files, @link_merge_files), 'pipelines with changed datasource ran ok - original and withdrawn files all exist';
my @new_files;
$element_count = 0;
foreach my $eid (map { $_->id } @{ get_elements($test_pipelinesetup_3->datasource) }) {
    $element_count++;
    my @output_subdirs = output_subdirs($eid, 3);
    my $step_index = 0;
    foreach my $suffix ('txt', 'txt.step_one') {
        push(@new_files, file(@output_subdirs, ($step_index + 1) . '_' . $expected_link_step_names[$step_index], "merged.$suffix"));
        $step_index++;
    }
}
is $element_count, 1, 'datasource 3 still has just 1 active element after withdrawing a parent element';
isnt $new_files[0], $link_merge_files[0], 'the merge pipeline has produced a new merge file since one of the parent elements got withdrawn';
ok handle_pipeline(@new_files), 'new merge files all exist after changing datasource';

my $element1 = VRPipe::DataElement->create(datasource => 2, result => { paths => [$base_files[2]->stringify] });
my $element2 = VRPipe::DataElement->create(datasource => 3, result => { paths => [map { $_->stringify } @base_files[1, 5, 9]], group => '50' });
is_deeply [$element1->withdrawn, $element2->withdrawn], [1, 1], 'elements correctly withdrawn after changing datasource';

# test that a child pipeline that deletes its input does not delete until the
# parent has finished with that file
my $output_dir_parent  = get_output_dir('parent_test_pipeline');
my $output_dir_deleter = get_output_dir('deleter_test_pipeline');

my $parent_pipeline = VRPipe::Pipeline->create(name => 'parent_pipeline', description => 'simple test pipeline with 3 steps, where the 3rd uses the output of the first');
VRPipe::StepAdaptor->create(pipeline => $parent_pipeline, to_step => 1, adaptor_hash => { one_input   => { data_element => 0 } });
VRPipe::StepAdaptor->create(pipeline => $parent_pipeline, to_step => 2, adaptor_hash => { delay_input => { one_output   => 1 } });
my $delay_step = VRPipe::Step->create(
    name               => 'delay_step',
    inputs_definition  => { delay_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'step two input file') },
    body_sub           => sub { my $self = shift; $self->dispatch_vrpipecode('sleep(60);', $self->new_requirements(memory => 500, time => 1)); },
    outputs_definition => {},
    post_process_sub => sub { return 1 },
    description      => 'the delay step'
);
VRPipe::StepAdaptor->create(pipeline => $parent_pipeline, to_step => 3, adaptor_hash => { two_input => { one_output => 1 } });
$parent_pipeline->add_step(VRPipe::Step->create(name => 'test_step_one'));
$parent_pipeline->add_step($delay_step);
$parent_pipeline->add_step(VRPipe::Step->create(name => 'test_step_two'));

my $parent_setup = VRPipe::PipelineSetup->create(
    name       => 'parent setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => $ds
    ),
    output_root => $output_dir_parent,
    pipeline    => $parent_pipeline,
    options     => {
        all_option => 'foo',
        one_option => 50
    }
);
my $delete_setup = VRPipe::PipelineSetup->create(
    name       => 'deleter setup',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'parent setup[1]',
        options => {}
    ),
    output_root => $output_dir_deleter,
    pipeline    => $link_pipeline,
    options     => {
        all_option    => 'bar',
        one_option    => 100,
        cleanup       => 1,
        delete_inputs => 1
    }
);
handle_pipeline();

my (@expected_files, @not_expected_files);
$element_id = 1;
foreach my $in ('file2.txt', 'file3.txt') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 4);
    push(@not_expected_files, file(@output_subdirs, '1_' . $expected_base_step_names[0], "$in.step_one"));
    push(@expected_files,     file(@output_subdirs, '3_' . $expected_base_step_names[1], "$in.step_one.step_two"));
}
foreach my $eid (map { $_->id } @{ get_elements($delete_setup->datasource) }) {
    my @output_subdirs = output_subdirs($eid, 5);
    push(@not_expected_files, file(@output_subdirs, '1_' . $expected_link_step_names[0], "merged.txt"));
    push(@expected_files,     file(@output_subdirs, '2_' . $expected_link_step_names[1], "merged.txt.step_one"));
}
is_deeply [scalar(@expected_files), scalar(@not_expected_files)], [4, 4], 'parent and deleter pipelines had the expected elements';
ok handle_pipeline(@expected_files), 'expected files exist after running a parent pipeline and child pipeline that deletes the parent output';
my $deleted = 0;
foreach my $path (@not_expected_files) {
    $deleted++ if (!-e $path);
}
# If vrpipe datasource is such that the child pipeline completes before parent
# step 2 completes, child will have deleted parent step 1 output before parent
# step 3 starts. Parent will then be forced to repeat step 1 to satisfy step 3's
# input needs. This is bad, and the next test would fail with only 2 deleted
# files (cleaned up files of child pipeline). If vrpipe datasource ensures that
# each parent element completes the pipeline before child starts, then no steps
# get repeated, and parent step 1 output remains deleted.
is $deleted, 4, 'the files that should have gotten deleted were deleted';

# # change metadata and detect changes
# VRPipe::File->create(path => $base_files[5]->stringify)->add_metadata({one_meta => '20'}, replace_data => 1);
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
# my $element3 = VRPipe::DataElement->create(datasource => 3, result => { paths => [ map { $_->stringify } @base_files[5,9] ], group => '50' });
# is $element3->withdrawn, 1, 'elements correctly withdrawn after metadata change';
#
# # if filter option one_meta#50 had not been set, a further data element
# # would have been created. check that we have the expected number of elements
# my $elements = get_elements($test_pipelinesetup_3->datasource);
# is scalar(@$elements), 1, 'group_by_metadata/filter option worked correctly - element not created';

finish;
