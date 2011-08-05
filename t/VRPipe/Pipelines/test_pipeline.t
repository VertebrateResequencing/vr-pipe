#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 8;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPipelines;
}

my $output_dir = get_output_dir('test_pipeline');
my $output_dir_clean = get_output_dir('test_pipeline_cleanup');

ok my $pipeline = VRPipe::Pipeline->get(name => 'test_pipeline'), 'able to get the test_pipeline pipeline';
my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(test_step_one test_step_two test_step_three test_step_four);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

$pipeline = VRPipe::Pipeline->get(name => 'test_pipeline');
@s_names = ();
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps after a second retrieval';

my $test_pipelinesetup = VRPipe::PipelineSetup->get(name => 'my test pipeline setup',
                                                    datasource => VRPipe::DataSource->get(type => 'fofn',
                                                                                          method => 'all',
                                                                                          source => file(qw(t data datasource.fofn2))),
                                                    output_root => $output_dir,
                                                    pipeline => $pipeline,
                                                    options => {all_option => 'foo',
                                                                one_option => 50,
                                                                four_option => 'bar',
                                                                cleanup => 0});

# also test with cleanup defaulting to true
my $test_pipelinesetup_clean = VRPipe::PipelineSetup->get(name => 'my test pipeline setup',
                                                          datasource => VRPipe::DataSource->get(type => 'fofn',
                                                                                                method => 'all',
                                                                                                source => file(qw(t data datasource.fofn2))),
                                                          output_root => $output_dir_clean,
                                                          pipeline => $pipeline,
                                                          options => {all_option => 'foo',
                                                                      one_option => 50,
                                                                      four_option => 'bar'});

my (@output_files, @final_files, @deleted_files);
my $element_id = 0;
foreach my $in ('file.txt', 'file2.txt', 'file3.txt') {
    $element_id++;
    my $step_index = 0;
    foreach my $suffix ('step_one', 'step_one.step_two', 'step_one.step_two.step_three', 'step_one.step_two.step_three.step_four') {
        push(@output_files, file($output_dir, output_subdirs($element_id), $expected_step_names[$step_index], "$in.$suffix"));
        
        if ($suffix eq 'step_one.step_two.step_three.step_four') {
            push(@final_files, file($output_dir_clean, output_subdirs($element_id), $expected_step_names[$step_index], "$in.$suffix"));
        }
        else {
            push(@deleted_files, file($output_dir_clean, output_subdirs($element_id), $expected_step_names[$step_index], "$in.$suffix"));
        }
        
        $step_index++;
    }
}
ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

my $ofile = VRPipe::File->get(path => file($output_dir, output_subdirs(3), 'test_step_four', 'file3.txt.step_one.step_two.step_three.step_four'));
my $ometa = $ofile->metadata;
my $o2meta = VRPipe::File->get(path => file($output_dir, output_subdirs(2), 'test_step_four', 'file2.txt.step_one.step_two.step_three.step_four'))->metadata;
is_deeply [$ometa->{one_meta}, $ometa->{two_meta}, $ometa->{three_meta}, $o2meta->{three_meta}, $ometa->{four_meta}], [50, 'body_decided_two_option', 'no_three_meta', 'StepOption_default_decided_three_option', 'bar'], 'metadata of one of the final output files was as expected';

my $existing_files = 0;
foreach my $file (@deleted_files) {
    $existing_files += -e $file ? 1 : 0;
}
is $existing_files, 0, 'all but the final files were deleted from the run with cleanup enabled';

my $expected_output = "3: a text file\n3: with two lines\n";
is_deeply [scalar($ofile->slurp), scalar(VRPipe::File->get(path => file($output_dir_clean, output_subdirs(3), 'test_step_four', 'file3.txt.step_one.step_two.step_three.step_four'))->slurp)], [$expected_output, $expected_output], 'both runs of the pipeline gave good output files';

done_testing;
exit;