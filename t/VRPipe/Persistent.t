#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Path::Class qw(file dir);
use File::Copy;

BEGIN {
    use Test::Most tests => 106;
    use VRPipeTest;
    
    use_ok('VRPipe::Persistent');
    use_ok('VRPipe::Persistent::Schema');
}

# some quick basic tests for all the core domain classes
my @schedulers;
ok $schedulers[0] = VRPipe::Scheduler->get(type => 'lsf'), 'created a Scheduler using get()';
is_deeply [$schedulers[0]->id, $schedulers[0]->type], [1, 'lsf'], 'scheduler1 has the expected fields';
ok $schedulers[1] = VRPipe::Scheduler->get(type => 'local'), 'created another Scheduler using get()';
is_deeply [$schedulers[1]->id, $schedulers[1]->type], [2, 'local'], 'scheduler2 has the expected fields';
ok my $default_type = VRPipe::Scheduler->default_type, 'could get a default type';
ok my $default_output_root = VRPipe::Scheduler->default_output_root, 'could get a default output root';
ok $schedulers[2] = VRPipe::Scheduler->get(), 'created another Scheduler using get() with no args';
is_deeply [$schedulers[2]->id, $schedulers[2]->type, $schedulers[2]->output_root], [1, $default_type, $default_output_root], 'scheduler3 has default fields';
$schedulers[2]->start_scheduler;

my $output_dir = dir($schedulers[2]->output_root, 'persistent_test_output');
$schedulers[2]->remove_tree($output_dir);
$schedulers[2]->make_path($output_dir);

my @files;
my $input1_path = file($output_dir, 'input1.txt');
open(my $fh, '>', $input1_path) or die "Could not write to $input1_path\n";
print $fh "input1_line1\ninput2_line2\n";
close($fh);
ok $files[0] = VRPipe::File->get(path => $input1_path, type => 'txt', metadata => {foo => 'bar'}), 'created a File using get()';
my $output1_path = file($output_dir, 'output1.txt');
ok $files[1] = VRPipe::File->get(path => $output1_path, type => 'txt'), 'created another File using get()';

my @ids;
ok $ids[0] = VRPipe::StepIODefinition->get(type => 'bam', description => 'step_1 bam input'), 'created a InputDefinition using get()';

my @steps;
ok $steps[0] = VRPipe::Step->get(name => 'step_1',
                                 inputs_definition => { static_input => $files[0],
                                                        dynamic_input => $ids[0] },
                                 body_sub => sub {
                                                    my $self = shift;
                                                    my $ofile = $self->output_file(output_key => 'step1_output', basename => 'output1.txt', type => 'txt');
                                                    my $fh = $ofile->openw();
                                                    print $fh "step1output\n";
                                                    $ofile->close();
                                                 },
                                 outputs_definition => { step1_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'step1_output file') },
                                 post_process_sub => sub { return 1 },
                                 description => 'the first step'), 'created a Step using get()';
is_deeply [$steps[0]->id, $steps[0]->description], [1, 'the first step'], 'step1 has the expected fields';
undef $steps[0];
ok $steps[0] = VRPipe::Step->get(name => 'step_1'), 'got a Step using get(name => )';
is_deeply [$steps[0]->id,
           $steps[0]->description,
           $steps[0]->inputs_definition->{static_input}->path], [1, 'the first step', $input1_path], 'step1 still has the expected fields';

ok my $first_step = VRPipe::Step->get(id => 1), 'step1 could be gotten by id';
is $first_step->description, 'the first step', 'it has the correct description';
is $first_step->description('the 1st step'), 'the first step', 'description does not seem like it was changed prior to an update';
is $steps[0]->description, 'the first step', 'indeed, other instances are not affected either';
$first_step->update;
is $steps[0]->description, 'the 1st step', 'all instances are affected after an update';

my @jobs;
my $epoch_time = time();
ok $jobs[0] = VRPipe::Job->get(cmd => 'echo "job1";'), 'created a Job using get() with no dir';
is_deeply [$jobs[0]->id, $jobs[0]->cmd, $jobs[0]->dir, $jobs[0]->running], [1, 'echo "job1";', cwd(), 0], 'job1 has the expected fields';
undef $jobs[0];
$jobs[0] = VRPipe::Job->get(cmd => 'echo "job1";');
cmp_ok $jobs[0]->creation_time->epoch, '>=', $epoch_time, 'creation time defaulted to just now';
is $jobs[0]->start_time, undef, 'start_time defaults to undef';
$jobs[1] = VRPipe::Job->get(cmd => 'echo "job2";');
$jobs[1]->pid(33);
$jobs[1]->user('foo');
$jobs[1]->host('local');
$jobs[1]->update;
undef $jobs[1];
$jobs[1] = VRPipe::Job->get(id => 2);
is_deeply [$jobs[1]->id, $jobs[1]->pid, $jobs[1]->user, $jobs[1]->host], [2, 33, 'foo', 'local'], 'you can set multiple fields followed by a single update';
$jobs[1]->pid(undef);
$jobs[1]->update;
undef $jobs[1];
$jobs[1] = VRPipe::Job->get(id => 2);
is $jobs[1]->pid, undef, 'you can set a field to undef';

my @reqs;
ok $reqs[0] = VRPipe::Requirements->get(memory => 2000, time => 6), 'created a Requirments using get() with only memory and time';
undef $reqs[0];
$reqs[0] = VRPipe::Requirements->get(memory => 2000, time => 6);
is_deeply [$reqs[0]->id, $reqs[0]->memory, $reqs[0]->time, $reqs[0]->cpus, $reqs[0]->tmp_space, $reqs[0]->local_space, $reqs[0]->custom], [1, 2000, 6, 1, 0, 0, {}], 'reqs1 has the expected fields';
ok $reqs[1] = VRPipe::Requirements->get(memory => 2000, time => 6, cpus => 2, tmp_space => 500, local_space => 0, custom => {loo => 'car'}), 'created a Requirments using fully specified get()';
undef $reqs[1];
$reqs[1] = VRPipe::Requirements->get(memory => 2000, time => 6, cpus => 2, tmp_space => 500, local_space => 0, custom => {loo => 'car'});
is_deeply [$reqs[1]->id, $reqs[1]->memory, $reqs[1]->time, $reqs[1]->cpus, $reqs[1]->tmp_space, $reqs[1]->local_space, $reqs[1]->custom->{loo}], [2, 2000, 6, 2, 500, 0, 'car'], 'reqs2 has the expected fields and is a seperate new entry in the db';

my @ds;
ok $ds[0] = VRPipe::DataSource->get(type => 'list', method => 'all', source => 't/data/datasource.fivelist'), 'created a DataSource using get()';
undef $ds[0];
$ds[0] = VRPipe::DataSource->get(type => 'list', method => 'all', source => 't/data/datasource.fivelist');
is_deeply [$ds[0]->id, $ds[0]->type, $ds[0]->method, $ds[0]->source], [1, 'list', 'all', 't/data/datasource.fivelist'], 'ds1 has the expected fields';
$ds[1] = VRPipe::DataSource->get(type => 'list', method => 'all', source => 't/data/datasource.onelist');
is_deeply [$ds[1]->id, $ds[1]->type, $ds[1]->method, $ds[1]->source], [2, 'list', 'all', 't/data/datasource.onelist'], 'ds2 has the expected fields';

my @de;
# create
$ds[0]->elements;
foreach my $de_num (1..5) {
    # get
    push(@de, VRPipe::DataElement->get(id => $de_num));
}
is_deeply [$de[2]->id, $de[2]->datasource->id, $de[2]->result->{line}], [3, 1, 'fed_result_3'], 'de3 has the expected fields';

my @pipelines;
ok $pipelines[0] = VRPipe::Pipeline->get(name => 'p1', description => 'first test pipeline'), 'created a Pipeline using get()';
undef $pipelines[0];
$pipelines[0] = VRPipe::Pipeline->get(name => 'p1', description => 'first test pipeline');
is_deeply [$pipelines[0]->id, $pipelines[0]->name, $pipelines[0]->description], [1, 'p1', 'first test pipeline'], 'pipeline1 has the expected fields';

# create some more steps we can chain together into a proper pipeline
VRPipe::Step->get(name => "step_2",
                  inputs_definition => { step2_input => VRPipe::StepIODefinition->get(type => 'txt', description => 'step_2 input') },
                  body_sub => sub {
                                    my $ofile = shift->output_file(output_key => 'step2_output', basename => 'step2_output.txt', type => 'txt');
                                    my $fh = $ofile->openw();
                                    print $fh "step2output\n";
                                    $ofile->close();
                                   },
                  outputs_definition => { step2_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'step2_output file') },
                  post_process_sub => sub { return 1 });
VRPipe::Step->get(name => "step_3",
                  inputs_definition => { step3_input => VRPipe::StepIODefinition->get(type => 'txt', description => 'step_3 input') },
                  body_sub => sub { 
                                    my $self = shift;
                                    my $ofile = $self->output_file(output_key => 'step3_output', basename => 'step3_output.txt', type => 'txt')->path;
                                    $self->dispatch(["sleep 5; echo step3output > $ofile", $self->new_requirements(memory => 50, time => 1)]);
                                    $self->dispatch(["sleep 4;", $self->new_requirements(memory => 50, time => 1)]);
                                    $self->dispatch(["sleep 3;", $self->new_requirements(memory => 50, time => 1)]);
                                   },
                  outputs_definition => { step3_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'step3_output file') },
                  post_process_sub => sub { return 1 });
VRPipe::Step->get(name => "step_4",
                  inputs_definition => { step4_input => VRPipe::StepIODefinition->get(type => 'txt', description => 'step_4 input') },
                  body_sub => sub {
                                    my $ofile = shift->output_file(output_key => 'step4_output', basename => 'step4_basename.txt', type => 'txt');
                                    my $fh = $ofile->openw();
                                    print $fh "step4output\n";
                                    $ofile->close();
                                   },
                  outputs_definition => { step4_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'step4_output file') },
                  post_process_sub => sub { return 1 });
VRPipe::Step->get(name => "step_5",
                  inputs_definition => { step5_input => VRPipe::StepIODefinition->get(type => 'txt', description => 'step_5 input') },
                  body_sub => sub {
                                    my $ofile = shift->output_file(output_key => 'step5_output', basename => 'step5_output.txt', type => 'txt');
                                    my $fh = $ofile->openw();
                                    print $fh "step5output\n";
                                    $ofile->close();
                                   },
                  outputs_definition => { step5_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'step5_output file') },
                  post_process_sub => sub { return 1 });
foreach my $step_num (2..5) {
    # get
    push(@steps, VRPipe::Step->get(id => $step_num));
}
is_deeply [$steps[2]->id, $steps[2]->inputs_definition->{step3_input}->type, $steps[2]->outputs_definition->{step3_output}->description], [3, 'txt', 'step3_output file'], 'step3 has the expected fields';

# create a 5 step pipeline by creating stepmembers
my @stepms;
my $step_num = 0;
foreach my $step (@steps) {
    # create
    VRPipe::StepMember->get(step => $step, pipeline => $pipelines[0], step_number => ++$step_num);
}
for (1..5) {
    # get
    push(@stepms, VRPipe::StepMember->get(id => $_));
}
is_deeply [$stepms[2]->id, $stepms[2]->step->id, $stepms[2]->pipeline->id, $stepms[2]->step_number], [3, 3, 1, 3], 'stepmember3 has the expected fields';

my @setups;
my $step1_bam_input = file($output_dir, 'input.bam');
copy(file(qw(t data file.bam)), $step1_bam_input) || die "Could not copy to $step1_bam_input\n";
ok $setups[0] = VRPipe::PipelineSetup->get(name => 'ps1', datasource => $ds[0], output_root => $output_dir, pipeline => $pipelines[0], options => {dynamic_input => "$step1_bam_input", baz => 'loman'}), 'created a PipelineSetup using get()';
undef $setups[0];
$setups[0] = VRPipe::PipelineSetup->get(id => 1);
is_deeply [$setups[0]->id, $setups[0]->datasource->id, $setups[0]->pipeline->id, $setups[0]->options->{baz}, $setups[0]->options->{dynamic_input}], [1, 1, 1, 'loman', $step1_bam_input], 'pipelinesetup1 has the expected fields when retrievied with just id';
undef $setups[0];
$setups[0] = VRPipe::PipelineSetup->get(name => 'ps1', datasource => $ds[0], output_root => $output_dir, pipeline => $pipelines[0], options => {dynamic_input => "$step1_bam_input", baz => 'loman'});
is_deeply [$setups[0]->id, $setups[0]->datasource->id, $setups[0]->pipeline->id, $setups[0]->options->{baz}], [1, 1, 1, 'loman'], 'pipelinesetup1 has the expected fields when retrieved with a full spec';


ok my $first_setup = VRPipe::PipelineSetup->get(id => 1), 'pipelinesetup1 could be gotten by id';
is $first_setup->datasource->id, 1, 'it has the correct datasource';
$first_setup->datasource($ds[1]); # *** though we shouldn't be able to change an is_key value
$first_setup->update;
is $first_setup->datasource->id, 2, 'datasource could be changed';
is $setups[0]->datasource->id, 2, 'and since we did an update, other instances are affected as well';
$first_setup->datasource($ds[0]);
$first_setup->update;

my @stepstates;
ok $stepstates[0] = VRPipe::StepState->get(stepmember => $stepms[0], dataelement => $de[0], pipelinesetup => $setups[0]), 'created a StepState using get()';
undef $stepstates[0];
$stepstates[0] = VRPipe::StepState->get(id => 1);
is_deeply [$stepstates[0]->id, $stepstates[0]->stepmember->id, $stepstates[0]->dataelement->id, $stepstates[0]->pipelinesetup->id, $stepstates[0]->complete], [1, 1, 1, 1, 0], 'stepstate1 has the expected fields';
$stepstates[1] = VRPipe::StepState->get(stepmember => $stepms[0], dataelement => $de[1], pipelinesetup => $setups[0]);

my @subs;
throws_ok { VRPipe::Submission->get(job => $jobs[0], stepstate => $stepstates[0]) } qr/Attribute \(requirements\) is required/, 'requirements is required when created a Submission, even though it is not a key';
ok $subs[0] = VRPipe::Submission->get(job => $jobs[0], stepstate => $stepstates[0], requirements => $reqs[0]), 'Submission could be made without required scheduler, which can default';
undef $subs[0];
$subs[0] = VRPipe::Submission->get(id => 1);
is_deeply [$subs[0]->id, $subs[0]->job->cmd, $subs[0]->stepstate->id, $subs[0]->requirements->memory, $subs[0]->retries, $subs[0]->scheduled, $subs[0]->done], [1, 'echo "job1";', 1, 2000, 0, undef, 0], 'submission1 has the expected fields';
ok my $defaulted_scheduler = $subs[0]->scheduler, 'submission1 has a scheduler generated by default';
is $defaulted_scheduler->id, 1, 'The scheduler was the default one';
is $subs[0]->memory, 2000, 'requirement methods pass through to the requriements object';
ok $subs[0]->extra_memory(500), 'could set some extra memory';
is $subs[0]->memory, 2500, 'now our requirements object has more memory';
is $subs[0]->requirements->id, 3, 'and it is actually a new row in the db';
$subs[0]->update;
is $subs[0]->id, 1, 'yet the submission id did not change';
$subs[0]->extra_memory();
is $subs[0]->memory, 3500, 'extra_memory defaults to adding 1000MB';
$subs[1] = VRPipe::Submission->get(job => $jobs[0], stepstate => $stepstates[0], requirements => $reqs[1]);
is $subs[1]->id, 1, 'a Submission created with only the reqs differing results in the same submission';
undef $subs[1];
$subs[1] = VRPipe::Submission->get(id => 1);
is $subs[1]->requirements->id, 2, 'but getting it with a certain requirements changed the requirements for this submission in the db';
$subs[1] = VRPipe::Submission->get(job => $jobs[0], stepstate => $stepstates[1], requirements => $reqs[1]);

throws_ok { VRPipe::PersistentArray->get() } qr/needs id or members/, 'get() for PersistentArray fails with no args';
throws_ok { VRPipe::PersistentArray->get(id => 1, members => \@subs) } qr/both id and members cannot be supplied/, 'get() for PersistentArray fails with both id and members supplied';
ok my $subs_array = VRPipe::PersistentArray->get(members => \@subs), 'created a PersistentArray using get(members => [...])';
is_deeply [$subs_array->id, $subs_array->members->[0]->id, $subs_array->members->[1]->id], [1, 1, 2], 'the created PArray has the correct contents';
undef $subs_array;
ok $subs_array = VRPipe::PersistentArray->get(id => 1), 'got a PersistentArray using get(id => 1)';
is_deeply [$subs_array->id, $subs_array->members->[0]->id, $subs_array->members->[1]->id], [1, 1, 2], 'the gotten PArray has the correct contents';
ok $subs_array = VRPipe::PersistentArray->get(members => \@subs), 'created a PersistentArray using the same set of members)';
is_deeply [$subs_array->id, $subs_array->members->[0]->id, $subs_array->members->[1]->id], [2, 1, 2], 'the created PArray has a new id, but the same contents otherwise';
is $subs_array->member(2)->id, 2, 'member() works given an index';
undef $subs_array;
$subs_array = VRPipe::PersistentArray->get(id => 1);
is $subs_array->member(2)->id, 2, 'member() works given an index when members() has not been called';

# steps can be created by requesting a name corresponding to a pre-written
# class in VRPipe::Steps::*
ok my $prewritten_step = VRPipe::Step->get(name => "md5_file_production"), 'able to get a pre-written step';

my %heartbeats;
# running jobs directly
my $tempdir = $jobs[1]->tempdir();
$jobs[2] = VRPipe::Job->get(cmd => qq[echo "job3"; sleep 3; perl -e 'print "foo\n"'], dir => $tempdir);
is $jobs[2]->heartbeat_interval(1), 1, 'heartbeat_interval of a job can be set directly and transiently';
$epoch_time = time();
$jobs[2]->run;
is_deeply [$jobs[2]->finished, $jobs[2]->running, $jobs[2]->ok, $jobs[2]->exit_code], [1, 0, 1, 0], 'test job status got updated correctly for an ok job';
ok $jobs[2]->pid, 'pid was set';
ok $jobs[2]->host, 'host was set';
ok $jobs[2]->user, 'user was set';
my $start_time = $jobs[2]->start_time->epoch;
my $end_time = $jobs[2]->end_time->epoch;
my $heartbeat_time = $jobs[2]->heartbeat->epoch;
my $ok = $start_time >= $epoch_time && $start_time <= $epoch_time + 1;
ok $ok, 'start_time is correct';
$ok = $end_time > $start_time && $end_time <= $start_time + 4;
ok $ok, 'end_time is correct';
$ok = $heartbeat_time > $start_time && $heartbeat_time <= $end_time;
ok $ok, 'time of last heartbeat correct';
ok my $stdout_file = $jobs[2]->stdout_file, 'got a stdout file';
is $stdout_file->slurp(chomp => 1), 'job3foo', 'stdout file had correct contents';
ok my $stderr_file = $jobs[2]->stderr_file, 'got a stderr file';
is $stderr_file->slurp(chomp => 1), '', 'stderr file was empty';

$jobs[2]->run;
is $jobs[2]->end_time->epoch, $end_time, 'running a job again does nothing';
ok $jobs[2]->reset_job, 'could reset a job';
is_deeply [$jobs[2]->finished, $jobs[2]->running, $jobs[2]->ok, $jobs[2]->exit_code, $jobs[2]->pid, $jobs[2]->host, $jobs[2]->user, $jobs[2]->heartbeat, $jobs[2]->start_time, $jobs[2]->end_time],
          [0, 0, 0, undef, undef, undef, undef, undef, undef, undef], 'after reset, job has cleared values';
my $child_pid = fork();
if ($child_pid) {
    sleep(1);
    my $cmd_pid = $jobs[2]->pid;
    kill(9, $cmd_pid);
    waitpid($child_pid, 0);
    is_deeply [$jobs[2]->finished, $jobs[2]->running, $jobs[2]->ok, $jobs[2]->exit_code], [1, 0, 0, 9], 'test job status got updated correctly for a job that was killed externally';
    is $jobs[2]->stdout_file->slurp(chomp => 1), 'job3', 'stdout file had correct contents';
    is $jobs[2]->stderr_file->slurp(chomp => 1), '', 'stderr file was empty';
}
else {
    $jobs[2]->run;
    exit(0);
}

$jobs[3] = VRPipe::Job->get(cmd => qq[echo "job4"; perl -e 'die "bar\n"'], dir => $tempdir);
$jobs[3]->heartbeat_interval(1);
$jobs[3]->run;
is_deeply [$jobs[3]->finished, $jobs[3]->running, $jobs[3]->ok, $jobs[3]->exit_code, $jobs[3]->heartbeat], [1, 0, 0, 65280, undef], 'test job status got updated correctly for a job that dies internally';
is $jobs[3]->stdout_file->slurp(chomp => 1), 'job4', 'stdout file had correct contents';
is $jobs[3]->stderr_file->slurp(chomp => 1), 'bar', 'stderr file had the correct contents';
throws_ok { $jobs[3]->run } qr/could not be run because it was not in the pending state/, 'run() on a failed job results in a throw';

# running jobs via the scheduler
$jobs[4] = VRPipe::Job->get(cmd => qq[perl -e 'foreach (1..5) { print "\$_\n"; sleep(1); }'], dir => $output_dir);
my $test_sub = VRPipe::Submission->get(job => $jobs[4], stepstate => $stepstates[0], requirements => $reqs[0]);
ok my $scheduled_id = $schedulers[2]->submit(submission => $test_sub), 'submit to the scheduler worked';
my $tlimit = time() + 1800;
wait_until_done($test_sub);
ok $test_sub->done, 'submission ran to completion';
is std_to_str($test_sub->job_stdout), "1\n2\n3\n4\n5\n", 'the submissions job did really run correctly';
$test_sub = VRPipe::Submission->get(job => $jobs[4], stepstate => $stepstates[1], requirements => $reqs[0]);
$schedulers[2]->submit(submission => $test_sub);
wait_until_done($test_sub);
is std_to_str($test_sub->job_stdout), "\n", 'running the same job in a different submission does not really rerun the job';

# job arrays
my @subs_array;
my @test_jobs;
for my $i (1..5) {
    my $output_dir = dir($schedulers[2]->output_root, 'test_output', $i);
    push(@test_jobs, VRPipe::Job->get(cmd => qq[perl -e 'foreach (1..9) { print "\$_\n"; sleep(1); } print \$\$, "\n"'], dir => $output_dir));
    push(@subs_array, VRPipe::Submission->get(job => $test_jobs[-1], stepstate => $stepstates[0], requirements => $reqs[0]));
}
ok $scheduled_id = $schedulers[2]->submit(array => \@subs_array), 'submit to the scheduler worked with an array';
throws_ok { $schedulers[2]->submit(array => \@subs_array); } qr/failed to claim all submissions/, 'trying to submit the same submissions again causes a throw';
%heartbeats = ();
wait_until_done(@subs_array);
my $good_outputs = 0;
foreach my $sub (@subs_array) {
    $good_outputs++ if std_to_str($sub->job_stdout) eq join("\n", 1..9)."\n".$sub->job->pid."\n";
}
is $good_outputs, scalar(@test_jobs), 'stdout files of all arrayed jobs had the correct contents';
my $good_beats = 0;
while (my ($sub_id, $hhash) = each %heartbeats) {
    my $beats = keys %{$hhash};
    my $num = keys %{$hhash};
    $good_beats++ if ($num >= 3 && $num <= 4);
}
is $good_beats, 5, 'each arrayed job had the correct number of heartbeats';

# stress testing
my ($t1, $l1);
$tlimit = time() + 1800;
SKIP: {
    skip "stress tests not enabled", 3 unless $ENV{VRPIPE_STRESSTESTS};
    
    my @subs_array;
    my @test_jobs;
    start_clock(__LINE__);
    for my $i (1..1000) {
        my $output_dir = dir($schedulers[2]->output_root, 'test_output', $i);
        push(@test_jobs, VRPipe::Job->get(cmd => qq[perl -e 'foreach (1..300) { print "\$_\n"; sleep(1); } print \$\$, "\n"'], dir => $output_dir));
        push(@subs_array, VRPipe::Submission->get(job => $test_jobs[-1], stepstate => $stepstates[0], requirements => $reqs[0]));
    }
    #                lustre nfs lustre
    lap(__LINE__); # 26 58 34
    ok my $scheduled_id = $schedulers[2]->submit(array => \@subs_array, heartbeat_interval => 30), 'submit to the scheduler worked with an array';
    lap(__LINE__); # 41 48 49
    %heartbeats = ();
    wait_until_done(@subs_array);
    lap(__LINE__); # 1046 2613 1818
    my $good_outputs = 0;
    foreach my $sub (@subs_array) {
        $good_outputs++ if $sub->job_stdout eq join("\n", 1..300)."\n".$sub->job->pid."\n";
    }
    lap(__LINE__); # 12 192 171
    is $good_outputs, scalar(@test_jobs), 'stdout files of all arrayed jobs had the correct contents';
    my $good_beats = 0;
    while (my ($sub_id, $hhash) = each %heartbeats) {
        my $beats = keys %{$hhash};
        $good_beats++ if ($beats >= 7 && $beats <= 10);
    }
    is $good_beats, 1000, 'each arrayed job had the correct number of heartbeats';
}

$schedulers[2]->stop_scheduler;
done_testing;
exit;

sub wait_until_done {
    my $loops = 0;
    my $all_done;
    while (1) {
        $all_done = 1;
        foreach my $sub (@_) {
            if (! $sub->job->finished) {
                $all_done = 0;
                my $heartbeat = $sub->job->heartbeat || next;
                $heartbeats{$sub->id}->{$heartbeat->epoch}++;
            }
        }
        last if $all_done;
        last if ++$loops > 1000;
        last if time() > $tlimit;
        sleep(1);
    }
    
    if ($all_done) {
        foreach my $sub (@_) {
            $sub->update_status;
        }
    }
}

sub start_clock {
    $l1 = shift;
    $t1 = time();
}

sub lap {
    my $l2 = shift;
    my $t2 = time();
    warn "Going from line $l1..$l2 took ", $t2 - $t1, " seconds\n";
    $t1 = time();
    $l1 = $l2 + 1;
}

sub std_to_str {
    my $pars = shift;
    $pars->next_record;
    return join("\n", @{$pars->parsed_record})."\n";
}