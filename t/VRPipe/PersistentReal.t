#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use File::Spec;

BEGIN {
    use Test::Most tests => 106;
    
    use_ok('VRPipe::Persistent');
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPersistentReal;
}

# some quick basic tests for all the core domain classes
my @steps;
ok $steps[0] = VRPipe::Step->get(name => 'coderef_1',
                                 inputs_sub => sub { return "inputs1" },
                                 body_sub => sub { return "body1" },
                                 post_process_sub => sub { return "post1" },
                                 outputs_sub => sub { return "outputs1" },
                                 description => 'the first step'), 'created a Step using get()';
is_deeply [$steps[0]->id, &{$steps[0]->body_sub}(), $steps[0]->description], [1, 'body1', 'the first step'], 'step1 has the expected fields';
undef $steps[0];
ok $steps[0] = VRPipe::Step->get(name => 'coderef_1'), 'got a Step using get(name => )';
is_deeply [$steps[0]->id, &{$steps[0]->body_sub}(), $steps[0]->description], [1, 'body1', 'the first step'], 'step1 still has the expected fields';

ok my $first_step = VRPipe::Step->get(id => 1), 'step1 could be gotten by id';
is $first_step->description, 'the first step', 'it has the correct description';
is $first_step->description('the 1st step'), 'the first step', 'description does not seem like it was changed prior to an update';
is $steps[0]->description, 'the first step', 'indeed, other instances are not affected either';
$first_step->update;
is $steps[0]->description, 'the 1st step', 'all instances are affected after an update';

my @schedulers;
ok $schedulers[0] = VRPipe::Scheduler->get(type => 'LSF', output_root => '/foo'), 'created a Scheduler using get()';
is_fields [qw/id type output_root/], $schedulers[0], [1, 'LSF', '/foo'], 'scheduler1 has the expected fields';
ok $schedulers[1] = VRPipe::Scheduler->get(type => 'local', output_root => '/bar'), 'created another Scheduler using get()';
is_fields [qw/id type output_root/], $schedulers[1], [2, 'local', '/bar'], 'scheduler2 has the expected fields';
ok my $default_type = VRPipe::Scheduler->default_type, 'could get a default type';
ok my $default_output_root = VRPipe::Scheduler->default_output_root, 'could get a default output root';
ok $schedulers[2] = VRPipe::Scheduler->get(), 'created another Scheduler using get() with no args';
is_fields [qw/id type output_root/], $schedulers[2], [3, $default_type, $default_output_root], 'scheduler3 has default fields';

my @jobs;
my $epoch_time = time();
ok $jobs[0] = VRPipe::Job->get(cmd => 'echo "job1";'), 'created a Job using get() with no dir';
is_fields [qw/id cmd dir running/], $jobs[0], [1, 'echo "job1";', cwd(), 0], 'job1 has the expected fields';
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
is_fields [qw/id pid user host/], $jobs[1], [2, 33, 'foo', 'local'], 'you can set multiple fields followed by a single update';
$jobs[1]->pid(undef);
$jobs[1]->update;
undef $jobs[1];
$jobs[1] = VRPipe::Job->get(id => 2);
is $jobs[1]->pid, undef, 'you can set a field to undef';

my @reqs;
ok $reqs[0] = VRPipe::Requirements->get(memory => 2000, time => 6), 'created a Requirments using get() with only memory and time';
undef $reqs[0];
$reqs[0] = VRPipe::Requirements->get(memory => 2000, time => 6);
is_fields [qw/id memory time cpus tmp_space local_space custom/], $reqs[0], [1, 2000, 6, 1, 0, 0, ''], 'reqs1 has the expected fields';
ok $reqs[1] = VRPipe::Requirements->get(memory => 2000, time => 6, cpus => 2, tmp_space => 500, local_space => 0, custom => 'foo'), 'created a Requirments using fully specified get()';
undef $reqs[1];
$reqs[1] = VRPipe::Requirements->get(memory => 2000, time => 6, cpus => 2, tmp_space => 500, local_space => 0, custom => 'foo');
is_fields [qw/id memory time cpus tmp_space local_space custom/], $reqs[1], [2, 2000, 6, 2, 500, 0, 'foo'], 'reqs2 has the expected fields and is a seperate new entry in the db';

my @ds;
ok $ds[0] = VRPipe::DataSource->get(module => 'VRPipe::DataSources::FOFN', method => 'all', source => 't/data/datasource.fofn'), 'created a DataSource using get()';
undef $ds[0];
$ds[0] = VRPipe::DataSource->get(module => 'VRPipe::DataSources::FOFN', method => 'all', source => 't/data/datasource.fofn');
is_fields [qw/id module method source/], $ds[0], [1, 'VRPipe::DataSources::FOFN', 'all', 't/data/datasource.fofn'], 'ds1 has the expected fields';
$ds[1] = VRPipe::DataSource->get(module => 'VRPipe::DataSources::FOFN', method => 'all', source => 't/data/datasource2.fofn');
is_fields [qw/id module method source/], $ds[1], [2, 'VRPipe::DataSources::FOFN', 'all', 't/data/datasource2.fofn'], 'ds2 has the expected fields';

my @de;
foreach my $de_num (1..5) {
    # create
    VRPipe::DataElement->get(datasource => $ds[0], result => "result_$de_num");
}
foreach my $de_num (1..5) {
    # get
    push(@de, VRPipe::DataElement->get(datasource => $ds[0], result => "result_$de_num"));
}
is_deeply [$de[2]->id, $de[2]->datasource->id, $de[2]->result], [3, 1, 'result_3'], 'de3 has the expected fields';

my @pipelines;
ok $pipelines[0] = VRPipe::Pipeline->get(name => 'p1', description => 'first test pipeline'), 'created a Pipeline using get()';
undef $pipelines[0];
$pipelines[0] = VRPipe::Pipeline->get(name => 'p1', description => 'first test pipeline');
is_fields [qw/id name description/], $pipelines[0], [1, 'p1', 'first test pipeline'], 'pipeline1 has the expected fields';

# create some more steps we can chain together into a proper pipeline
VRPipe::Step->get(name => "coderef_2",
                  inputs_sub => sub { return "inputs2" },
                  body_sub => sub { return "body2" },
                  post_process_sub => sub { return "post2" },
                  outputs_sub => sub { return "outputs2" });
VRPipe::Step->get(name => "coderef_3",
                  inputs_sub => sub { return "inputs3" },
                  body_sub => sub { return "body3" },
                  post_process_sub => sub { return "post3" },
                  outputs_sub => sub { return "outputs3" });
VRPipe::Step->get(name => "coderef_4",
                  inputs_sub => sub { return "inputs4" },
                  body_sub => sub { return "body4" },
                  post_process_sub => sub { return "post4" },
                  outputs_sub => sub { return "outputs4" });
VRPipe::Step->get(name => "coderef_5",
                  inputs_sub => sub { return "inputs5" },
                  body_sub => sub { return "body5" },
                  post_process_sub => sub { return "post5" },
                  outputs_sub => sub { return "outputs5" });
foreach my $step_num (2..5) {
    # get
    push(@steps, VRPipe::Step->get(id => $step_num));
}
is_deeply [$steps[2]->id, &{$steps[2]->body_sub}()], [3, 'body3'], 'step3 has the expected fields';

# create a 5 step pipeline by creating stepmembers
my @stepms;
foreach my $step (@steps) {
    # create
    VRPipe::StepMember->get(step => $step, pipeline => $pipelines[0]);
}
for (1..5) {
    # get
    push(@stepms, VRPipe::StepMember->get(id => $_));
}
is_deeply [$stepms[2]->id, $stepms[2]->step->id, $stepms[2]->pipeline->id], [3, 3, 1], 'stepmember3 has the expected fields';

my @setups;
ok $setups[0] = VRPipe::PipelineSetup->get(name => 'ps1', datasource => $ds[0], output_root => '/foo/bar', pipeline => $pipelines[0], options => {foo => 'bar', baz => 'loman'}), 'created a PipelineSetup using get()';
undef $setups[0];
$setups[0] = VRPipe::PipelineSetup->get(id => 1);
is_deeply [$setups[0]->id, $setups[0]->datasource->id, $setups[0]->pipeline->id, $setups[0]->options->{baz}], [1, 1, 1, 'loman'], 'pipelinesetup1 has the expected fields when retrievied with just id';
undef $setups[0];
$setups[0] = VRPipe::PipelineSetup->get(name => 'ps1', datasource => $ds[0], output_root => '/foo/bar', pipeline => $pipelines[0], options => {foo => 'bar', baz => 'loman'});
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
is $defaulted_scheduler->id, 3, 'The scheduler was the default one';
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
ok $jobs[2]->reset, 'could reset a job';
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
my %heartbeats;
my $output_dir = File::Spec->catdir($schedulers[2]->output_root, 'test_output');
$jobs[4] = VRPipe::Job->get(cmd => qq[perl -e 'foreach (1..5) { print "\$_\n"; sleep(1); }'], dir => $output_dir);
my $test_sub = VRPipe::Submission->get(job => $jobs[4], stepstate => $stepstates[0], requirements => $reqs[0]);
ok my $scheduled_id = $schedulers[2]->submit(submission => $test_sub), 'submit to the scheduler worked';
wait_until_done($test_sub);
ok $test_sub->done, 'submission ran to completion';
my $job4_stdout_file = $jobs[4]->stdout_file;
is $job4_stdout_file->slurp(chomp => 1), '12345', 'the submissions job did really run correctly';
$test_sub = VRPipe::Submission->get(job => $jobs[4], stepstate => $stepstates[1], requirements => $reqs[0]);
$schedulers[2]->submit(submission => $test_sub);
wait_until_done($test_sub);
is $jobs[4]->stdout_file, $job4_stdout_file, 'running the same job in a different submission does not really rerun the job';

# job arrays
my @subs_array;
my @test_jobs;
for my $i (1..5) {
    my $output_dir = File::Spec->catdir($schedulers[2]->output_root, 'test_output', $i);
    push(@test_jobs, VRPipe::Job->get(cmd => qq[perl -e 'foreach (1..9) { print "\$_\n"; sleep(1); } print \$\$, "\n"'], dir => $output_dir));
    push(@subs_array, VRPipe::Submission->get(job => $test_jobs[-1], stepstate => $stepstates[0], requirements => $reqs[0]));
}
ok $scheduled_id = $schedulers[2]->submit(array => \@subs_array), 'submit to the scheduler worked with an array';
throws_ok { $schedulers[2]->submit(array => \@subs_array); } qr/failed to claim all submissions/, 'trying to submit the same submissions again causes a throw';
%heartbeats = ();
wait_until_done(@subs_array);
my $good_outputs = 0;
foreach my $job (@test_jobs) {
    $good_outputs++ if $job->stdout_file->slurp(chomp => 1) eq join('', 1..9).$job->pid;
}
is $good_outputs, scalar(@test_jobs), 'stdout files of all arrayed jobs had the correct contents';
my $good_beats = 0;
while (my ($sub_id, $hhash) = each %heartbeats) {
    my $beats = keys %{$hhash};
    $good_beats++ if keys %{$hhash} == 3;
}
is $good_beats, 5, 'each arrayed job had the correct number of heartbeats';

# manager
ok my $manager = VRPipe::Manager->get(), 'can get a manager with no args';
is $manager->id, 1, 'manager always has an id of 1';
$pipelines[1] = VRPipe::Pipeline->get(name => 'p2', description => 'second pipeline');
$setups[1] = VRPipe::PipelineSetup->get(name => 'ps2', datasource => $ds[0], output_root => '/flar/blar', pipeline => $pipelines[1]);
my @manager_setups = $manager->setups;
is @manager_setups, 2, 'setups() returns the correct number of PipelineSetups';
@manager_setups = $manager->setups(pipeline_name => 'p2');
is @manager_setups, 1, 'setups() returns the correct number of PipelineSetups when pipeline_name supplied';
$manager->trigger;

# stress testing
my ($t1, $l1);
SKIP: {
    skip "stress tests not enabled", 3 unless $ENV{VRPIPE_STRESSTESTS};
    
    my @subs_array;
    my @test_jobs;
    start_clock(__LINE__);
    for my $i (1..1000) {
        my $output_dir = File::Spec->catdir($schedulers[2]->output_root, 'test_output', $i);
        push(@test_jobs, VRPipe::Job->get(cmd => qq[perl -e 'foreach (1..300) { print "\$_\n"; sleep(1); } print \$\$, "\n"'], dir => $output_dir));
        push(@subs_array, VRPipe::Submission->get(job => $test_jobs[-1], stepstate => $stepstates[0], requirements => $reqs[0]));
    }
    lap(__LINE__); # 26
    ok $scheduled_id = $schedulers[2]->submit(array => \@subs_array, heartbeat_interval => 30), 'submit to the scheduler worked with an array';
    lap(__LINE__); # 41
    %heartbeats = ();
    wait_until_done(@subs_array);
    lap(__LINE__); # 1046
    my $good_outputs = 0;
    foreach my $job (@test_jobs) {
        $good_outputs++ if $job->stdout_file->slurp(chomp => 1) eq join('', 1..300).$job->pid;
    }
    lap(__LINE__); # 12
    is $good_outputs, scalar(@test_jobs), 'stdout files of all arrayed jobs had the correct contents';
    my $good_beats = 0;
    while (my ($sub_id, $hhash) = each %heartbeats) {
        my $beats = keys %{$hhash};
        $good_beats++ if keys %{$hhash} == 10;
    }
    is $good_beats, 1000, 'each arrayed job had the correct number of heartbeats';
}

done_testing;
exit;

sub wait_until_done {
    my $loops = 0;
    while (1) {
        my $all_done = 1;
        foreach my $sub (@_) {
            if (! $sub->done) {
                $all_done = 0;
                #last;
                my $heartbeat = $sub->job->heartbeat || next;
                $heartbeats{$sub->id}->{$heartbeat->epoch}++;
            }
        }
        last if $all_done;
        last if ++$loops > 1000;
        sleep(1);
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