#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;

BEGIN {
    use Test::Most tests => 26;
    
    use_ok('VRPipe::Persistent');
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPersistentReal;
}

# some quick basic tests for all the core domain classes
my @steps;
ok $steps[0] = VRPipe::Step->get(code => 'step1code', description => 'the first step'), 'created a Step using get()';
is_fields [qw/id code description/], $steps[0], [1, 'step1code', 'the first step'], 'step1 has the expected fields';
undef $steps[0];
ok $steps[0] = VRPipe::Step->get(code => 'step1code', description => 'the first step'), 'got a Step using get()';
is_fields [qw/id code description/], $steps[0], [1, 'step1code', 'the first step'], 'step1 still has the expected fields';

my @schedulers;
ok $schedulers[0] = VRPipe::Scheduler->get(module => 'VRPipe::Schedulers::LSF'), 'created a Scheduler using get()';
is_fields [qw/id module/], $schedulers[0], [1, 'VRPipe::Schedulers::LSF'], 'scheduler1 has the expected fields';
ok $schedulers[1] = VRPipe::Scheduler->get(module => 'VRPipe::Schedulers::Local'), 'created another Scheduler using get()';
is_fields [qw/id module/], $schedulers[1], [2, 'VRPipe::Schedulers::Local'], 'scheduler2 has the expected fields';

my @jobs;
my $epoch_time = time();
ok $jobs[0] = VRPipe::Job->get(cmd => 'echo "job1";'), 'created a Job using get() with no dir';
is_fields [qw/id cmd dir running/], $jobs[0], [1, 'echo "job1";', cwd(), 0], 'job1 has the expected fields';
undef $jobs[0];
$jobs[0] = VRPipe::Job->get(cmd => 'echo "job1";');
cmp_ok $jobs[0]->creation_time->epoch, '>=', $epoch_time, 'creation time defaulted to just now';
is $jobs[0]->start_time, undef, 'start_time defaults to undef';

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
    VRPipe::DataElement->get(datasource => $ds[1], result => "result_$de_num");
}
foreach my $de_num (1..5) {
    # get
    push(@de, VRPipe::DataElement->get(datasource => $ds[1], result => "result_$de_num"));
}
is_deeply [$de[2]->id, $de[2]->datasource->id, $de[2]->result], [3, 2, 'result_3'], 'de3 has the expected fields';

my @pipelines;
ok $pipelines[0] = VRPipe::Pipeline->get(name => 'p1', description => 'first test pipeline'), 'created a Pipeline using get()';
undef $pipelines[0];
$pipelines[0] = VRPipe::Pipeline->get(name => 'p1', description => 'first test pipeline');
is_fields [qw/id name description/], $pipelines[0], [1, 'p1', 'first test pipeline'], 'pipeline1 has the expected fields';

foreach my $step_num (2..5) {
    # create
    VRPipe::Step->get(code => "coderef_$step_num");
}
foreach my $step_num (2..5) {
    # get
    push(@steps, VRPipe::Step->get(code => "coderef_$step_num"));
}
is_fields [qw/id code/], $steps[2], [3, 'coderef_3'], 'step3 has the expected fields';

my @stepms;
foreach my $step (@steps) {
    # create
    VRPipe::StepMember->get(step => $step, pipeline => $pipelines[0]);
}
foreach my $step (@steps) {
    # get
    push(@stepms, VRPipe::StepMember->get(step => $step, pipeline => $pipelines[0]));
}
is_deeply [$stepms[2]->id, $stepms[2]->step->id, $stepms[2]->pipeline->id], [3, 3, 1], 'stepmember3 has the expected fields';


#my @subs;
#ok $subs[0] = VRPipe::Submission->get(), 'created a Submission using get()';


done_testing;
exit;