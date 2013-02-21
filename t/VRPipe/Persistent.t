#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Path::Class qw(file dir);
use File::Copy;
use Parallel::ForkManager;
use DateTime;
use EV;
use AnyEvent;

BEGIN {
    use Test::Most tests => 151;
    use VRPipeTest;
    
    use_ok('VRPipe::Persistent');
    use_ok('VRPipe::Persistent::Schema');
}

# some quick basic tests for all the core domain classes
my @schedulers;
ok my $default_type = VRPipe::Scheduler->default_type, 'could get a default type';
ok $schedulers[0] = VRPipe::Scheduler->create(type => 'lsf'), 'created a Scheduler using create()';
is $schedulers[0]->type, 'lsf', 'scheduler1 has the expected type';
ok $schedulers[1] = VRPipe::Scheduler->create(type => 'local'), 'created another Scheduler using create()';
is $schedulers[1]->type, 'local', 'scheduler2 has the expected type';
ok my $default_output_root = VRPipe::Scheduler->default_output_root, 'could get a default output root';
ok $schedulers[2] = VRPipe::Scheduler->get(), 'created another Scheduler using get() with no args';
is_deeply [$schedulers[2]->id, $schedulers[2]->type, $schedulers[2]->output_root], [$default_type eq 'lsf' ? $schedulers[0]->id : $schedulers[1]->id, $default_type, $default_output_root], 'scheduler3 has default fields';
$schedulers[2]->start_scheduler;

my $output_dir = dir($schedulers[2]->output_root, 'persistent_test_output');
$schedulers[2]->remove_tree($output_dir);
$schedulers[2]->make_path($output_dir);

my @files;
my $input1_path = file($output_dir, 'input1.txt');
open(my $fh, '>', $input1_path) or die "Could not write to $input1_path\n";
print $fh "input1_line1\ninput2_line2\n";
close($fh);
ok $files[0] = VRPipe::File->create(path => $input1_path, type => 'txt', metadata => { foo => 'bar' }), 'created a File using create()';
my $output1_path = file($output_dir, 'output1.txt');
ok $files[1] = VRPipe::File->create(path => $output1_path, type => 'txt'), 'created another File using create()';

$files[0]->disconnect;
my $post_disconnect = VRPipe::File->get(id => 1);
my $pd_worked = $post_disconnect && $post_disconnect->path() eq $input1_path;
ok $pd_worked, 'we automatically reconnected to the database after a disconnect';

my @ids;
ok $ids[0] = VRPipe::StepIODefinition->create(type => 'bam', description => 'step_1 bam input'), 'created a InputDefinition using create()';

my @steps;
ok $steps[0] = VRPipe::Step->create(
    name              => 'step_1',
    inputs_definition => {
        static_input  => $files[0],
        dynamic_input => $ids[0]
    },
    body_sub => sub {
        my $self  = shift;
        my $ofile = $self->output_file(output_key => 'step1_output', basename => 'output1.txt', type => 'txt');
        my $fh    = $ofile->openw();
        print $fh "step1output\n";
        $ofile->close();
    },
    outputs_definition => { step1_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step1_output file') },
    post_process_sub   => sub            { return 1 },
    description        => 'the first step'
  ),
  'created a Step using create()';
is_deeply [$steps[0]->id, $steps[0]->description], [1, 'the first step'], 'step1 has the expected fields';
undef $steps[0];
ok $steps[0] = VRPipe::Step->get(name => 'step_1'), 'got a Step using get(name => )';
is_deeply [$steps[0]->id, $steps[0]->description, $steps[0]->inputs_definition->{static_input}->path], [1, 'the first step', $input1_path], 'step1 still has the expected fields';

ok my $first_step = VRPipe::Step->get(id => 1), 'step1 could be gotten by id';
is $first_step->description, 'the first step', 'it has the correct description';
is $first_step->description('the 1st step'), 'the first step', 'description does not seem like it was changed prior to an update';
is $steps[0]->description, 'the first step', 'indeed, other instances are not affected either';
$first_step->update;
$steps[0]->reselect_values_from_db;
is $steps[0]->description, 'the 1st step', 'a different instance is affected after an update and reselect';

my @jobs;
my $epoch_time = time();
ok $jobs[0] = VRPipe::Job->create(cmd => 'echo "job1";'), 'created a Job using create() with no dir';
is_deeply [$jobs[0]->id, $jobs[0]->cmd, $jobs[0]->dir], [1, 'echo "job1";', cwd()], 'job1 has the expected fields';
undef $jobs[0];
$jobs[0] = VRPipe::Job->get(cmd => 'echo "job1";');
is $jobs[0]->start_time, undef, 'start_time defaults to undef';
$jobs[1] = VRPipe::Job->create(cmd => 'echo "job2";');
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
ok $reqs[0] = VRPipe::Requirements->create(memory => 2000, time => 6), 'created a Requirments using create() with only memory and time';
undef $reqs[0];
$reqs[0] = VRPipe::Requirements->get(memory => 2000, time => 6);
is_deeply [$reqs[0]->id, $reqs[0]->memory, $reqs[0]->time, $reqs[0]->cpus, $reqs[0]->tmp_space, $reqs[0]->local_space, $reqs[0]->custom], [1, 2000, 21600, 1, 0, 0, {}], 'reqs1 has the expected fields';
ok $reqs[1] = VRPipe::Requirements->create(memory => 2000, time => 6, cpus => 2, tmp_space => 500, local_space => 0, custom => { loo => 'car' }), 'created a Requirments using fully specified create()';
undef $reqs[1];
$reqs[1] = VRPipe::Requirements->get(memory => 2000, time => 6, cpus => 2, tmp_space => 500, local_space => 0, custom => { loo => 'car' });
is_deeply [$reqs[1]->id, $reqs[1]->memory, $reqs[1]->time, $reqs[1]->cpus, $reqs[1]->tmp_space, $reqs[1]->local_space, $reqs[1]->custom->{loo}], [2, 2000, 21600, 2, 500, 0, 'car'], 'reqs2 has the expected fields and is a seperate new entry in the db';

my @ds;
ok $ds[0] = VRPipe::DataSource->create(type => 'list', method => 'all', source => 't/data/datasource.fivelist'), 'created a DataSource using create()';
undef $ds[0];
$ds[0] = VRPipe::DataSource->get(type => 'list', method => 'all', source => 't/data/datasource.fivelist');
is_deeply [$ds[0]->id, $ds[0]->type, $ds[0]->method, $ds[0]->source], [1, 'list', 'all', 't/data/datasource.fivelist'], 'ds1 has the expected fields';
$ds[1] = VRPipe::DataSource->create(type => 'list', method => 'all', source => 't/data/datasource.onelist');
is_deeply [$ds[1]->id, $ds[1]->type, $ds[1]->method, $ds[1]->source], [2, 'list', 'all', 't/data/datasource.onelist'], 'ds2 has the expected fields';

my @de;
# create
$ds[0]->elements;
foreach my $de_num (1 .. 5) {
    # get
    push(@de, VRPipe::DataElement->get(id => $de_num));
}
is_deeply [$de[2]->id, $de[2]->datasource->id, $de[2]->result->{line}], [3, 1, 'fed_result_3'], 'de3 has the expected fields';

my @pipelines;
ok $pipelines[0] = VRPipe::Pipeline->create(name => 'p1', description => 'first test pipeline'), 'created a Pipeline using create()';
undef $pipelines[0];
$pipelines[0] = VRPipe::Pipeline->get(name => 'p1', description => 'first test pipeline');
is_deeply [$pipelines[0]->id, $pipelines[0]->name, $pipelines[0]->description], [1, 'p1', 'first test pipeline'], 'pipeline1 has the expected fields';

# create some more steps we can chain together into a proper pipeline
VRPipe::Step->create(
    name              => "step_2",
    inputs_definition => { step2_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'step_2 input') },
    body_sub          => sub {
        my $ofile = shift->output_file(output_key => 'step2_output', basename => 'step2_output.txt', type => 'txt');
        my $fh = $ofile->openw();
        print $fh "step2output\n";
        $ofile->close();
    },
    outputs_definition => { step2_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step2_output file') },
    post_process_sub   => sub            { return 1 }
);
VRPipe::Step->create(
    name              => "step_3",
    inputs_definition => { step3_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'step_3 input') },
    body_sub          => sub {
        my $self = shift;
        my $ofile = $self->output_file(output_key => 'step3_output', basename => 'step3_output.txt', type => 'txt')->path;
        $self->dispatch(["sleep 5; echo step3output > $ofile", $self->new_requirements(memory => 50, time => 1)]);
        $self->dispatch(["sleep 4;",                           $self->new_requirements(memory => 50, time => 1)]);
        $self->dispatch(["sleep 3;",                           $self->new_requirements(memory => 50, time => 1)]);
    },
    outputs_definition => { step3_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step3_output file') },
    post_process_sub   => sub            { return 1 }
);
VRPipe::Step->create(
    name              => "step_4",
    inputs_definition => { step4_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'step_4 input') },
    body_sub          => sub {
        my $ofile = shift->output_file(output_key => 'step4_output', basename => 'step4_basename.txt', type => 'txt');
        my $fh = $ofile->openw();
        print $fh "step4output\n";
        $ofile->close();
    },
    outputs_definition => { step4_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step4_output file') },
    post_process_sub   => sub            { return 1 }
);
VRPipe::Step->create(
    name              => "step_5",
    inputs_definition => { step5_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'step_5 input') },
    body_sub          => sub {
        my $ofile = shift->output_file(output_key => 'step5_output', basename => 'step5_output.txt', type => 'txt');
        my $fh = $ofile->openw();
        print $fh "step5output\n";
        $ofile->close();
    },
    outputs_definition => { step5_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step5_output file') },
    post_process_sub   => sub            { return 1 }
);
foreach my $step_num (2 .. 5) {
    # get
    push(@steps, VRPipe::Step->get(id => $step_num));
}
is_deeply [$steps[2]->id, $steps[2]->inputs_definition->{step3_input}->type, $steps[2]->outputs_definition->{step3_output}->description], [3, 'txt', 'step3_output file'], 'step3 has the expected fields';

# create a 5 step pipeline by creating stepmembers
my @stepms;
my $step_num = 0;
foreach my $step (@steps) {
    # create
    VRPipe::StepMember->create(step => $step, pipeline => $pipelines[0], step_number => ++$step_num);
}
for (1 .. 5) {
    # get
    push(@stepms, VRPipe::StepMember->get(id => $_));
}
is_deeply [$stepms[2]->id, $stepms[2]->step->id, $stepms[2]->pipeline->id, $stepms[2]->step_number], [3, 3, 1, 3], 'stepmember3 has the expected fields';

is $pipelines[0]->num_steps, 5, 'num_steps worked';

my @setups;
my $step1_bam_input = file($output_dir, 'input.bam');
copy(file(qw(t data file.bam)), $step1_bam_input) || die "Could not copy to $step1_bam_input\n";
ok $setups[0] = VRPipe::PipelineSetup->create(name => 'ps1', datasource => $ds[0], output_root => $output_dir, pipeline => $pipelines[0], options => { dynamic_input => "$step1_bam_input", baz => 'loman' }), 'created a PipelineSetup using create()';
undef $setups[0];
$setups[0] = VRPipe::PipelineSetup->get(id => 1);
is_deeply [$setups[0]->id, $setups[0]->datasource->id, $setups[0]->pipeline->id, $setups[0]->options->{baz}, $setups[0]->options->{dynamic_input}], [1, 1, 1, 'loman', $step1_bam_input], 'pipelinesetup1 has the expected fields when retrievied with just id';
undef $setups[0];
$setups[0] = VRPipe::PipelineSetup->get(name => 'ps1', datasource => $ds[0], output_root => $output_dir, pipeline => $pipelines[0], options => { dynamic_input => "$step1_bam_input", baz => 'loman' });
is_deeply [$setups[0]->id, $setups[0]->datasource->id, $setups[0]->pipeline->id, $setups[0]->options->{baz}], [1, 1, 1, 'loman'], 'pipelinesetup1 has the expected fields when retrieved with a full spec';
# (we test desired_farm specifically because we have an around alteration of it)
is $setups[0]->desired_farm, undef, 'desired_farm defaults to undef';
$setups[0]->desired_farm('foo');
$setups[0]->update;
is $setups[0]->desired_farm, 'foo', 'desired_farm could be set and retrieved';

ok my $first_setup = VRPipe::PipelineSetup->get(id => 1), 'pipelinesetup1 could be gotten by id';
is $first_setup->datasource->id, 1, 'it has the correct datasource';
$first_setup->datasource($ds[1]); # *** though we shouldn't be able to change an is_key value
$first_setup->update;
is $first_setup->datasource->id, 2, 'datasource could be changed';
$setups[0]->reselect_values_from_db;
is $setups[0]->datasource->id, 2, 'and since we did an update, other instances are affected as well';
$first_setup->datasource($ds[0]);
$first_setup->update;

ok my $log = VRPipe::PipelineSetupLog->create(ps_id => 1, message => "foo bar", date => DateTime->now()), 'able to call PipelineSetupLog->create';
is $log->ps_id,   1,         'ps_id() worked';
is $log->message, 'foo bar', 'message() worked';
ok !$log->stack, 'strack() defaults to undef';
$first_setup->log_event('created and changed');
is_deeply [map { $_->message } $first_setup->logs()], ['foo bar', 'created and changed'], 'logs() returned all the logs for the first setup';
is_deeply [map { $_->message } $first_setup->logs(like => 'created%')], ['created and changed'], q[logs with like returned all the logs for the first setup with matching string];
sleep(2);
$first_setup->log_event('created and changed');
is_deeply [map { $_->message } $first_setup->logs()], ['foo bar', 'created and changed'], 'log_event with the same message does not work twice in a row';
$first_setup->log_event('car');
$first_setup->log_event('created and changed');
is_deeply [map { $_->message } $first_setup->logs()], ['foo bar', 'created and changed', 'car', 'created and changed'], 'log_event does allow the same message to be recorded more than once';
is_deeply [$log->de_id, $log->ss_id, $log->sub_id, $log->job_id], [0, 0, 0, 0], 'de_id, ss_id, sub_id, job_id default to 0';
$log = $first_setup->log_event("funion", record_stack => 1, dataelement => 1, stepstate => 2, submission => 3, job => 4);
is_deeply [$log->de_id, $log->ss_id, $log->sub_id, $log->job_id], [1, 2, 3, 4], 'de_id, ss_id, sub_id, job_id can all be set during log_event';
ok $log->stack, 'stack() can be set during log_event';
$first_setup->log_event("funion2", dataelement => 1);
$first_setup->log_event("funion3", dataelement => 2);
is_deeply [map { $_->message } $first_setup->logs(dataelement => 1)], ['funion', 'funion2'], 'logs() limited to a given dataelement returned only matching logs';
is_deeply [map { $_->message } $first_setup->logs(dataelement => 1, stepstate => 2)], ['funion'], 'logs() limited to a given dataelement and stepstate returned more restricted matching log';
is_deeply [map { $_->message } $first_setup->logs(dataelement => 1, stepstate => 2, include_undefined => 1)], ['foo bar', 'created and changed', 'car', 'created and changed', 'funion', 'funion2'], 'logs() limited to a given dataelement and stepstate with include_undefined on returned all logs except for de 2';

my @stepstates;
ok $stepstates[0] = VRPipe::StepState->create(stepmember => $stepms[0], dataelement => $de[0], pipelinesetup => $setups[0]), 'created a StepState using create()';
undef $stepstates[0];
$stepstates[0] = VRPipe::StepState->get(id => 1);
is_deeply [$stepstates[0]->id, $stepstates[0]->stepmember->id, $stepstates[0]->dataelement->id, $stepstates[0]->pipelinesetup->id, $stepstates[0]->complete], [1, 1, 1, 1, 0], 'stepstate1 has the expected fields';
$stepstates[1] = VRPipe::StepState->create(stepmember => $stepms[0], dataelement => $de[1], pipelinesetup => $setups[0]);

my @subs;
throws_ok { VRPipe::Submission->create(job => $jobs[0], stepstate => $stepstates[0]) } qr/Attribute \(requirements\) is required/, 'requirements is required when created a Submission, even though it is not a key';
ok $subs[0] = VRPipe::Submission->create(job => $jobs[0], stepstate => $stepstates[0], requirements => $reqs[0]), 'Submission could be made without required scheduler, which can default';
undef $subs[0];
$subs[0] = VRPipe::Submission->get(id => 1);
is_deeply [$subs[0]->id, $subs[0]->job->cmd, $subs[0]->stepstate->id, $subs[0]->requirements->memory, $subs[0]->retries, $subs[0]->_claim, $subs[0]->done], [1, 'echo "job1";', 1, 2000, 0, 0, 0], 'submission1 has the expected fields';
ok my $defaulted_scheduler = $subs[0]->scheduler, 'submission1 has a scheduler generated by default';
is $defaulted_scheduler->id, $default_type eq 'lsf' ? $schedulers[0]->id : $schedulers[1]->id, 'The scheduler was the default one';
is $subs[0]->memory, 2000, 'requirement methods pass through to the requirements object';
ok $subs[0]->extra_memory(500), 'could set some extra memory';
is $subs[0]->memory, 2500, 'now our requirements object has more memory';
is $subs[0]->requirements->id, 3, 'and it is actually a new row in the db';
$subs[0]->update;
is $subs[0]->id, 1, 'yet the submission id did not change';
$subs[0]->extra_memory();
is $subs[0]->memory, 3500, 'extra_memory defaults to adding 1000MB';
$subs[1] = VRPipe::Submission->create(job => $jobs[0], stepstate => $stepstates[0], requirements => $reqs[1]);
is $subs[1]->id, 1, 'a Submission created with only the reqs differing results in the same submission';
undef $subs[1];
$subs[1] = VRPipe::Submission->get(id => 1);
is $subs[1]->requirements->id, 2, 'but getting it with a certain requirements changed the requirements for this submission in the db';
$subs[1] = VRPipe::Submission->create(job => $jobs[0], stepstate => $stepstates[1], requirements => $reqs[1]);

throws_ok { VRPipe::PersistentArray->get() } qr/requires id or members/, 'get() for PersistentArray fails with no args';
throws_ok { VRPipe::PersistentArray->create() } qr/Validation failed/,   'create() for PersistentArray fails with no args';
throws_ok { VRPipe::PersistentArray->get(id => 1, members => \@subs) } qr/cannot supply both id and members/, 'get() for PersistentArray fails with both id and members supplied';
ok my $subs_array = VRPipe::PersistentArray->create(members => \@subs), 'created a PersistentArray using create(members => [...])';
is_deeply [$subs_array->id, ($subs_array->members)[0]->id, ($subs_array->members)[1]->id, $subs_array->member(1)->id, $subs_array->member(2)->id], [1, 1, 2, 1, 2], 'the created PArray has the correct contents';
undef $subs_array;
ok $subs_array = VRPipe::PersistentArray->get(id => 1), 'got a PersistentArray using get(id => 1)';
is_deeply [$subs_array->id, ($subs_array->members)[0]->id, ($subs_array->members)[1]->id], [1, 1, 2], 'the gotten PArray has the correct contents';
ok $subs_array = VRPipe::PersistentArray->create(members => \@subs), 'created a PersistentArray using the same set of members)';
is_deeply [$subs_array->id, ($subs_array->members)[0]->id, ($subs_array->members)[1]->id, $subs_array->member(1)->id, $subs_array->member(2)->id], [2, 3, 4, 1, 2], 'the created PArray has a new id, new persistentarray members, but the same contents otherwise';
ok $subs_array = VRPipe::PersistentArray->get(members => \@subs), 'got a PersistentArray using a set of members';
is $subs_array->id, 1, 'it had the id of the first matching PArray';
ok $subs_array = VRPipe::PersistentArray->get(members => [$subs[1], $subs[0]], any_order => 1), 'got a PersistentArray using get()';
is $subs_array->id, 1, 'using get() with members in a previously unused order still returns an existing PArray when any_order is true';
ok $subs_array = VRPipe::PersistentArray->get(members => [$subs[1], $subs[0]]), 'created a PersistentArray using get()';
is $subs_array->id, 3, 'using get() with members in a previously unused order returns a new PArray';
is_deeply [map { $_->id => ref($_) } $subs_array->member_instances], [2 => 'VRPipe::Submission', 1 => 'VRPipe::Submission'], 'member_instances() works';

# now that we have some submissions and pipelinesetup, make some stepstats and
# test that the multi-row get methods work
my @expected_stepstats = ([5, 1, 1, 1], [10, 2, 2, 2]);
my @stepstat_cols      = qw(memory time id submission);
my $search_args        = { step => 1 };
my $sub_id             = 1;
foreach my $es (@expected_stepstats) {
    VRPipe::StepStats->create(step => 1, pipelinesetup => 1, submission => $sub_id++, memory => $es->[0], time => $es->[1]);
}

my @got_stepstats = ();
my $rs            = Schema->resultset('StepStats')->search($search_args);
while (my $stepstat = $rs->next) {
    push(@got_stepstats, [$stepstat->memory, $stepstat->time, $stepstat->id, $stepstat->submission->id]);
}
is_deeply \@got_stepstats, \@expected_stepstats, 'using a manual resultset search got expected stepstat values';

@got_stepstats = ();
$rs            = VRPipe::StepStats->search_rs($search_args);
while (my $stepstat = $rs->next) {
    push(@got_stepstats, [$stepstat->memory, $stepstat->time, $stepstat->id, $stepstat->submission->id]);
}
is_deeply \@got_stepstats, \@expected_stepstats, 'using VRPipe::StepStats->search_rs and an $rs->next loop got expected stepstat values';

@got_stepstats = map { [$_->memory, $_->time, $_->id, $_->submission->id] } VRPipe::StepStats->search($search_args);
is_deeply \@got_stepstats, \@expected_stepstats, 'using VRPipe::StepStats->search in list context got expected stepstat values without an $rs->next loop';

is scalar(VRPipe::StepStats->search($search_args)), 2, 'using VRPipe::StepStats->search in scalar context gave the correct count of results';
is_deeply [map { $_->memory } VRPipe::StepStats->search($search_args, { rows => 1 })], [5], 'VRPipe::StepStats->search with a rows attribute works';

my $pager = VRPipe::StepStats->search_paged($search_args, {}, 1);
@got_stepstats = ();
my $pages = 0;
while (my $stepstats = $pager->next) {
    $pages++;
    foreach my $stepstat (@$stepstats) {
        push(@got_stepstats, [$stepstat->memory, $stepstat->time, $stepstat->id, $stepstat->submission->id]);
    }
}
is_deeply [\@got_stepstats, $pages], [\@expected_stepstats, 2], 'using VRPipe::StepStats->search_paged and an $pager->next loop got expected stepstat values';

@got_stepstats = @{ VRPipe::StepStats->get_column_values(\@stepstat_cols, $search_args) };
is_deeply \@got_stepstats, \@expected_stepstats, 'using get_column_values with multiple columns got expected stepstat values';

$pager         = VRPipe::StepStats->get_column_values_paged(\@stepstat_cols, $search_args, {}, 1);
@got_stepstats = ();
$pages         = 0;
while (my $vals = $pager->next) {
    $pages++;
    push(@got_stepstats, @$vals);
}
is_deeply [\@got_stepstats, $pages], [\@expected_stepstats, 2], 'using VRPipe::StepStats->get_column_values_paged and a $pager->next loop got expected stepstat values';

@got_stepstats = VRPipe::StepStats->get_column_values('memory', { step => VRPipe::Step->get(id => 1) });
is_deeply \@got_stepstats, [5, 10], 'using get_column_values with a single column and an instance in the search args got expected stepstat values';
is_deeply [VRPipe::StepStats->get_column_values('memory', { step => 1 }, { rows => 1 })], [5], 'get_column_values with a rows attribute works';

$pager         = VRPipe::StepStats->get_column_values_paged('memory', $search_args, {}, 1);
@got_stepstats = ();
$pages         = 0;
while (my $vals = $pager->next) {
    $pages++;
    push(@got_stepstats, $vals);
}
is_deeply [\@got_stepstats, $pages], [[[5], [10]], 2], 'using VRPipe::StepStats->get_column_values_paged with a single column and a $pager->next loop got expected stepstat values';

# test that we don't repeat pages when new rows are added between next() calls
$stepstates[2] = VRPipe::StepState->create(stepmember => $stepms[0], dataelement => $de[2], pipelinesetup => $setups[0]);
$subs[2] = VRPipe::Submission->create(job => $jobs[0], stepstate => $stepstates[2], requirements => $reqs[1]);

$pager         = VRPipe::StepStats->search_paged($search_args, {}, 1);
@got_stepstats = ();
$pages         = 0;
while (my $stepstats = $pager->next) {
    $pages++;
    foreach my $stepstat (@$stepstats) {
        push(@got_stepstats, [$stepstat->memory, $stepstat->time, $stepstat->id, $stepstat->submission->id]);
    }
    if ($pages == 1) {
        VRPipe::StepStats->create(step => 1, pipelinesetup => 1, submission => 3, memory => 15, time => 3);
    }
}
is_deeply [\@got_stepstats, $pages], [\@expected_stepstats, 2], 'search_paged and an $pager->next loop got expected (old) stepstat values, even when a new row got added between pages';

push(@expected_stepstats, [15, 3, 3, 3]);
$pager         = VRPipe::StepStats->search_paged($search_args, {}, 1);
@got_stepstats = ();
$pages         = 0;
while (my $stepstats = $pager->next) {
    $pages++;
    foreach my $stepstat (@$stepstats) {
        push(@got_stepstats, [$stepstat->memory, $stepstat->time, $stepstat->id, $stepstat->submission->id]);
    }
}
is_deeply [\@got_stepstats, $pages], [\@expected_stepstats, 3], 'a new pager gives us the new row we added';

$pager         = VRPipe::StepStats->search_paged({ %{$search_args}, 'me.id' => 1 }, {}, 1);
@got_stepstats = ();
$pages         = 0;
while (my $stepstats = $pager->next) {
    $pages++;
    foreach my $stepstat (@$stepstats) {
        push(@got_stepstats, [$stepstat->memory, $stepstat->time, $stepstat->id, $stepstat->submission->id]);
    }
}
is_deeply [\@got_stepstats, $pages], [[$expected_stepstats[0]], 1], 'search_paged auto-id_limiting does not stop user using their own id search terms';

# test that this works with something like Step, which has refs for vals
$pager = VRPipe::Step->get_column_values_paged(['body_sub', 'id'], { id => { '<=' => 2 } }, {}, 1);
my @got_steps = ();
$pages = 0;
while (my $vals = $pager->next) {
    $pages++;
    push(@got_steps, @$vals);
}
is_deeply [ref($got_steps[0]->[0]), $got_steps[0]->[1], scalar(@got_steps), $pages], ['CODE', 1, 2, 2], 'using VRPipe::Step->get_column_values_paged and a $pager->next loop got expected step values, where the refs were deflated';

$pager     = VRPipe::Step->get_column_values_paged('body_sub', { id => { '<=' => 2 } }, {}, 1);
@got_steps = ();
$pages     = 0;
while (my $vals = $pager->next) {
    $pages++;
    push(@got_steps, $vals);
}
is_deeply [ref($got_steps[0]->[0]), scalar(@got_steps), $pages], ['CODE', 2, 2], 'using VRPipe::StepStats->get_column_values_paged with a single column and a $pager->next loop got expected step values, where the refs were deflated';

use_ok('VRPipe::StepStatsUtil');
my $ssu = VRPipe::StepStatsUtil->new(step => VRPipe::Step->get(id => 1));
my @ssumm = $ssu->mean_memory;
is $ssumm[1], 10, 'StepStatsUtil mean_memory, which is implemented with get_column_values_paged, worked fine';
is_deeply [$ssu->percentile_memory(percent => 90, pipelinesetup => VRPipe::PipelineSetup->get(id => 1))], [3, 15], 'StepStatsUtil percentile_memory, which is implemented with get_column_values, worked fine';

my @job_args;
foreach my $i (1 .. 1000) {
    push(@job_args, { cmd => "fake_job $i", dir => '/fake_dir' });
}
VRPipe::Job->bulk_create_or_update(@job_args);
my $j_count = VRPipe::Job->search({ dir => '/fake_dir', exit_code => undef });
is $j_count, 1000, 'bulk_create_or_update worked when creating';
@job_args = ();
foreach my $i (1 .. 1000) {
    push(@job_args, { cmd => "fake_job $i", dir => '/fake_dir', exit_code => 0 });
}
VRPipe::Job->bulk_create_or_update(@job_args);
$j_count = VRPipe::Job->search({ dir => '/fake_dir', exit_code => undef });
my $j_count_exited = VRPipe::Job->search({ dir => '/fake_dir', exit_code => 0 });
is_deeply [$j_count, $j_count_exited], [0, 1000], 'bulk_create_or_update worked when updating, and no duplicate rows were created';
@job_args = ();
foreach my $i (1 .. 10) {
    push(@job_args, { cmd => "fake_job with running set $i", dir => '/fake_dir', start_time => DateTime->now, heartbeat => DateTime->now });
}
VRPipe::Job->bulk_create_or_update({ cmd => "fake_job withouth running set", dir => '/fake_dir' }, @job_args);
$j_count = VRPipe::Job->search({ heartbeat => { '!=' => undef } });
is $j_count, 10, 'bulk_create_or_update was able to create and set non keys when it was first supplied some args that lacked the non key column';

my $fm = Parallel::ForkManager->new(2);
for (1 .. 2) {
    $fm->start and next;
    
    @job_args = ();
    foreach my $i (1001 .. 2000) {
        push(@job_args, { cmd => "fake_job $i", dir => '/fake_dir_parallel' });
    }
    VRPipe::Job->bulk_create_or_update(@job_args);
    
    $fm->finish(0);
}
$fm->wait_all_children;
$j_count = VRPipe::Job->search({ dir => '/fake_dir_parallel' });
is $j_count, 1000, 'bulk_create_or_update worked when the same creation was requested by 2 different processes at the same time';

# steps can be created by requesting a name corresponding to a pre-written
# class in VRPipe::Steps::*
ok my $prewritten_step = VRPipe::Step->get(name => "md5_file_production"), 'able to get a pre-written step using get() instead of create()';

# so can pipelines; test that they get constructed ok
ok my $prewritten_pipeline = VRPipe::Pipeline->create(name => 'test_pipeline'), 'able to get a pre-written pipeline using create()';
is_deeply [map { $_->name } @{ $prewritten_pipeline->steps }], [qw(test_step_one test_step_two test_step_three test_step_four)], 'steps method worked';
is_deeply [map { $_->id } @{ $prewritten_pipeline->adaptors }], [1, 2, 3, 4], 'adaptors method worked';
is_deeply [map { $_->id } @{ $prewritten_pipeline->behaviours }], [1, 2, 3], 'behaviours method worked';

my %heartbeats;
# running jobs directly
my $tempdir = $jobs[1]->tempdir();
$jobs[2] = VRPipe::Job->create(cmd => qq[echo "job3"; sleep 3; perl -e 'print "foo\n"'], dir => $tempdir);
is $jobs[2]->heartbeat_interval(1), 1, 'heartbeat_interval of a job can be set directly and transiently';
$epoch_time = time();
$jobs[2]->run;
run_job($jobs[2]);
is_deeply [defined($jobs[2]->end_time), $jobs[2]->ok, $jobs[2]->exit_code], [1, 1, 0], 'test job status got updated correctly for an ok job';
ok $jobs[2]->pid,  'pid was set';
ok $jobs[2]->host, 'host was set';
ok $jobs[2]->user, 'user was set';
my $start_time     = $jobs[2]->start_time->epoch;
my $end_time       = $jobs[2]->end_time->epoch;
my $heartbeat_time = $jobs[2]->heartbeat->epoch;
my $ok             = $start_time >= $epoch_time && $start_time <= $epoch_time + 1;
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
run_job($jobs[2]);
is $jobs[2]->end_time->epoch, $end_time, 'running a job again does nothing';
ok $jobs[2]->reset_job, 'could reset a job';
is_deeply [$jobs[2]->ok, $jobs[2]->exit_code, $jobs[2]->pid, $jobs[2]->host, $jobs[2]->user, $jobs[2]->heartbeat, $jobs[2]->start_time, $jobs[2]->end_time], [0, undef, undef, undef, undef, undef, undef, undef], 'after reset, job has cleared values';
my $own_pid   = $$;
my $child_pid = fork();
if ($child_pid) {
    sleep(1);
    $jobs[2]->reselect_values_from_db;
    my $cmd_pid = $jobs[2]->pid;
    kill(9, $cmd_pid);
    waitpid($child_pid, 0);
    $jobs[2]->reselect_values_from_db;
    is_deeply [defined($jobs[2]->end_time), $jobs[2]->ok, $jobs[2]->exit_code], [1, 0, 9], 'test job status got updated correctly for a job that was killed externally';
    is $jobs[2]->stdout_file->slurp(chomp => 1), 'job3', 'stdout file had correct contents';
    is $jobs[2]->stderr_file->slurp(chomp => 1), '',     'stderr file was empty';
}
else {
    $jobs[2]->run;
    run_job($jobs[2]);
    exit(0);
}

$jobs[3] = VRPipe::Job->create(cmd => qq[echo "job4"; perl -e 'die "bar\n"'], dir => $tempdir);
$jobs[3]->heartbeat_interval(1);
$jobs[3]->run;
run_job($jobs[3]);
is_deeply [defined($jobs[3]->end_time), $jobs[3]->ok, $jobs[3]->exit_code =~ /65280|512/], [1, 0, 1], 'test job status got updated correctly for a job that dies internally';
is $jobs[3]->stdout_file->slurp(chomp => 1), 'job4', 'stdout file had correct contents';
is $jobs[3]->stderr_file->slurp(chomp => 1), 'bar',  'stderr file had the correct contents';
ok !$jobs[3]->run, 'run() on a failed job does not work';

# running jobs via a submission
my $t_before = time();
$jobs[4] = VRPipe::Job->create(cmd => qq[perl -e 'foreach (1..5) { print "\$_ ", time(), "\n"; sleep(1); }'], dir => $output_dir);
my $test_sub = VRPipe::Submission->create(job => $jobs[4], stepstate => $stepstates[0], requirements => $reqs[0]);
wait_until_done($test_sub);
my $t_after = time();
ok $test_sub->done, 'submission ran to completion';
my $std_to_str_first_time = std_to_str($test_sub->job_stdout);
like $std_to_str_first_time, qr/^1 (\d+)\n2 (\d+)\n3 (\d+)\n4 (\d+)\n5 (\d+)\n$/, 'the submissions job did really run correctly';
my ($job_time) = $std_to_str_first_time =~ /^1 (\d+)/;
my $job_time_correct = ($job_time >= $t_before) && ($job_time <= $t_after);
ok $job_time_correct, 'the job output we checked really was just made seconds ago';
$test_sub = VRPipe::Submission->create(job => $jobs[4], stepstate => $stepstates[1], requirements => $reqs[0]);
wait_until_done($test_sub);
is std_to_str($test_sub->job_stdout), "\n", 'running the same job in a different submission does not really rerun the job';

# (these tests were originally supposed to test running multiple submissions
#  simultaneously on the farm via the scheduler, but now they run them
#  sequentially on the local cpu... should fix this...)
my @subs_array;
my @test_jobs;
for my $i (1 .. 5) {
    my $output_dir = dir($schedulers[2]->output_root, 'test_output', $i);
    push(@test_jobs, VRPipe::Job->create(cmd => qq[perl -e 'foreach (1..9) { print "\$_\n"; sleep(1); } print \$\$, "\n"'], dir => $output_dir));
    push(@subs_array, VRPipe::Submission->create(job => $test_jobs[-1], stepstate => $stepstates[0], requirements => $reqs[0]));
}
%heartbeats = ();
wait_until_done(@subs_array);
my $good_outputs = 0;
foreach my $sub (@subs_array) {
    $good_outputs++ if std_to_str($sub->job_stdout) eq join("\n", 1 .. 9) . "\n" . $sub->job->pid . "\n";
}
is $good_outputs, scalar(@test_jobs), 'stdout files of all arrayed jobs had the correct contents';
my $good_beats = 0;
while (my ($sub_id, $hhash) = each %heartbeats) {
    my $beats = keys %{$hhash};
    my $num   = keys %{$hhash};
    $good_beats++ if ($num >= 3 && $num <= 4);
}
is $good_beats, 5, 'each arrayed job had the correct number of heartbeats';

$schedulers[2]->stop_scheduler;
done_testing;
exit;

sub run_job {
    my $job = shift;
    
    my $watcher = EV::timer 0, 2, sub {
        if ($job->end_time) {
            EV::unloop;
        }
    };
    
    EV::run;
}

sub wait_until_done {
    my @subs = @_;
    
    my $count = 0;
    my $sub;
    my $watcher = EV::timer 0, 2, sub {
        $count++;
        if ($count > 500) {
            EV::unloop;
        }
        
        if (!$sub || $sub->done || $sub->failed) {
            $sub = pop(@subs);
            unless ($sub) {
                EV::unloop;
                return;
            }
            $sub->claim_and_run;
        }
        
        my $job = $sub->job;
        if ($job->end_time) {
            $sub->update_status;
            $sub->release;
        }
        else {
            my $heartbeat = $sub->job->heartbeat || return;
            $heartbeats{ $sub->id }->{ $heartbeat->epoch }++;
        }
    };
    
    EV::run;
}

sub std_to_str {
    my $pars = shift;
    $pars->next_record;
    return join("\n", @{ $pars->parsed_record }) . "\n";
}
