#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 28;
    use VRPipeTest;
    use TestPipelines;
    
    use_ok('VRPipe::StepStatsUtil');
}

my $output_dir = get_output_dir('test_pipeline');
my $output_dir_clean = get_output_dir('test_pipeline_cleanup');

ok my $pipeline = VRPipe::Pipeline->create(name => 'test_pipeline'), 'able to get the test_pipeline pipeline';
my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(test_step_one test_step_two test_step_three test_step_four);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

$pipeline = VRPipe::Pipeline->create(name => 'test_pipeline');
@s_names = ();
my @s_nums = ();
my $last_step;
foreach my $stepmember ($pipeline->steps) {
    $last_step = $stepmember->step;
    push(@s_names, $last_step->name);
    push(@s_nums, $stepmember->step_number);
}
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps after a second retrieval';
is_deeply \@s_nums, [1, 2, 3, 4], 'the pipeline has the correct step numbers after a second retrieval';

my $test_pipelinesetup = VRPipe::PipelineSetup->create(name => 'my test pipeline setup',
                                                    datasource => VRPipe::DataSource->create(type => 'fofn',
                                                                                          method => 'all',
                                                                                          source => file(qw(t data datasource.fofn2))),
                                                    output_root => $output_dir,
                                                    pipeline => $pipeline,
                                                    options => {all_option => 'foo',
                                                                one_option => 50,
                                                                four_option => 'bar',
                                                                cleanup => 0});

# also test with cleanup defaulting to true
my $test_pipelinesetup_clean = VRPipe::PipelineSetup->create(name => 'my test pipeline setup',
                                                          datasource => VRPipe::DataSource->create(type => 'fofn',
                                                                                                method => 'all',
                                                                                                source => file(qw(t data datasource.fofn2))),
                                                          output_root => $output_dir_clean,
                                                          pipeline => $pipeline,
                                                          options => {all_option => 'foo',
                                                                      one_option => 50,
                                                                      four_option => 'bar'});

# before the step has ever been run, there is no recommended memory or time
ok my $ssu = VRPipe::StepStatsUtil->new(step => $last_step), 'able to make a StepStatsUtil object';
is $ssu->recommended_memory, undef, 'recommended_memory returns undef to start with';
is $ssu->recommended_time, undef, 'recommended_time returns undef to start with';

my (@output_files, @final_files, @deleted_files);
my $element_id = 0;
foreach my $in ('file.txt', 'file2.txt', 'file3.txt') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    my @clean_output_subdirs = output_subdirs($element_id, 2);
    
    my $step_index = 0;
    foreach my $suffix ('step_one', 'step_one.step_two', 'step_one.step_two.step_three', 'step_one.step_two.step_three.step_four') {
        push(@output_files, file(@output_subdirs, ($step_index+1).'_'.$expected_step_names[$step_index], "$in.$suffix"));
        
        if ($suffix eq 'step_one.step_two.step_three.step_four') {
            push(@final_files, file(@clean_output_subdirs, ($step_index+1).'_'.$expected_step_names[$step_index], "$in.$suffix"));
        }
        else {
            push(@deleted_files, file(@clean_output_subdirs, ($step_index+1).'_'.$expected_step_names[$step_index], "$in.$suffix"));
        }
        
        $step_index++;
    }
}
ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

my $ofile = VRPipe::File->create(path => file(output_subdirs(3), '4_test_step_four', 'file3.txt.step_one.step_two.step_three.step_four'));
my $ometa = $ofile->metadata;
my $o2meta = VRPipe::File->create(path => file(output_subdirs(2), '4_test_step_four', 'file2.txt.step_one.step_two.step_three.step_four'))->metadata;
is_deeply [$ometa->{one_meta}, $ometa->{two_meta}, $ometa->{three_meta}, $o2meta->{three_meta}, $ometa->{four_meta}], [50, 'body_decided_two_option', undef, 'StepOption_default_decided_three_option', 'bar'], 'metadata of one of the final output files was as expected';

my $existing_files = 0;
foreach my $file (@deleted_files) {
    $existing_files += -e $file ? 1 : 0;
}
is $existing_files, 0, 'all but the final files were deleted from the run with cleanup enabled';

my $expected_output = "3: a text file\n3: with two lines\n";
is_deeply [scalar($ofile->slurp), scalar(VRPipe::File->create(path => file(output_subdirs(3, 2), '4_test_step_four', 'file3.txt.step_one.step_two.step_three.step_four'))->slurp)], [$expected_output, $expected_output], 'both runs of the pipeline gave good output files';


# let's reset the final step of the cleaned-up pipeline and rerun, to test if
# we trigger a cascade of step resets when it discovers previous step output
# files have been deleted
my $last_stepstate = VRPipe::StepState->create(pipelinesetup => 2, dataelement => 3, stepmember => 4);
$last_stepstate->start_over;

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output after we did start_over on the last cleaned-up stepstate';


# let's also test manually starting a particular dataelement over from scratch
my $de_to_restart = VRPipe::DataElement->create(id => 3);
$de_to_restart->start_from_scratch($test_pipelinesetup_clean);
is_deeply [-e $final_files[-1], -e $final_files[-2]], [undef, 1], 'start_from_scratch deleted the final output file for the particular datalement, but not others';
ok handle_pipeline(@output_files, @final_files), 'pipeline ran and recreated files for the scratched dataelement';


# check that we can access job and scheduler std output/err files
my $submission = VRPipe::Submission->create(id => 1);
my $scheduler = VRPipe::Scheduler->create();
my $subm_dir = dir($scheduler->output_root, qw(b f a 2 VRPipe__Job__1));
is $submission->scheduler_stdout_file->path, file($subm_dir, 'scheduler_output_file'), 'scheduler_stdout_file was correct';
is $submission->scheduler_stderr_file->path, file($subm_dir, 'scheduler_error_file'), 'scheduler_stderr_file was correct';
is $submission->job_stdout_file->path, file($subm_dir, 'job_stdout'), 'job_stdout_file was correct';
is $submission->job_stderr_file->path, file($subm_dir, 'job_stderr'), 'job_stderr_file was correct';
my $pars = $submission->scheduler_stderr;
$pars->next_record;
is_deeply [@{$pars->parsed_record}], [], 'scheduler_stderr had no content';
my $parser = $submission->scheduler_stdout;
ok $parser->does('VRPipe::ParserRole'), 'scheduler_stdout returns a parser';


# let's test moving a mid-step output file and starting over a final step of the
# non-cleaned pipeline, to confirm that it does not redo the mid-step but uses
# the moved file
my $orig_file = VRPipe::File->create(path => $output_files[10]);
my $moved_file = VRPipe::File->create(path => $output_files[10].'.moved');
$orig_file->move($moved_file);
is_deeply [-e $output_files[10], -e $output_files[10].'.moved'], [undef, 1], 'we were able to move a step 3 output file';
$output_files[10] = $output_files[10].'.moved';
my %mtimes;
foreach my $i (0..$#output_files) {
    my $file = VRPipe::File->create(path => $output_files[$i]);
    $mtimes{$output_files[$i]} = $file->e ? $file->mtime : 0;
}
$last_stepstate = VRPipe::StepState->create(pipelinesetup => 1, dataelement => 3, stepmember => 4);
$last_stepstate->start_over;
my $orig_final = VRPipe::File->create(path => $output_files[11]);
my $new_final = $output_files[11];
$new_final =~ s/step_three/step_three.moved/;
$output_files[11] = $new_final;
$new_final = VRPipe::File->create(path => $new_final);

ok handle_pipeline(@output_files), 'pipeline ran and created all expected output after we did start_over on the last un-cleaned stepstate';
my $oks = 0;
foreach my $i (0..$#output_files) {
    my $file = VRPipe::File->create(path => $output_files[$i]);
    my $mtime = $file->e ? $file->mtime : 0;
    my $orig_mtime = $mtimes{$output_files[$i]} || 0;
    if ($i == 11) {
        $oks++ if $mtime ne $orig_mtime;
    }
    else {
        $oks++ if $mtime eq $orig_mtime;
    }
}
$orig_final->reselect_values_from_db;
$new_final->reselect_values_from_db;
is_deeply [$oks, -e $orig_file->path, $orig_final->e, $new_final->e], [12, undef, 0, 1], 'only the final step file we deleted was recreated (with a new name) - not the step 3 file we moved';

# now that the step has run a number of times, we should have recommended reqs
is $ssu->recommended_memory, 100, 'recommended_memory returns a new value at the end';
is $ssu->recommended_time, 1, 'recommended_time returns new value at the end';
is $ssu->recommended_time(pipelinesetup => $test_pipelinesetup_clean), 1, 'recommended_time can return values specific to a particular pipelinesetup';
is_deeply [($ssu->percentile_seconds(percent => 95))[0], ($ssu->percentile_seconds(percent => 95, pipelinesetup => $test_pipelinesetup_clean))[0]], [6, 3], 'percentile_seconds gives the correct counts the recommendations are based on';

finish;