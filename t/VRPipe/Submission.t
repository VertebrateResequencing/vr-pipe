#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class qw(file dir);
use Parallel::ForkManager;

BEGIN {
    use Test::Most tests => 8;
    use VRPipeTest;
}

# create some submissions (and everything needed to create them)
my $scheduler = VRPipe::Scheduler->create;
my $req       = VRPipe::Requirements->create(memory => 500, time => 3600);
my $pipeline  = VRPipe::Pipeline->create(name => 'test_pipeline');
my $ds        = VRPipe::DataSource->create(type => 'list', method => 'all', source => file(qw(t data datasource.onelist)));
my $ps        = VRPipe::PipelineSetup->create(name => 'ps', datasource => $ds, output_root => dir(qw(tmp)), pipeline => $pipeline, options => {});
$ps->active(0); # stop the server trying to do something with this ps
$ps->update;
$ds->elements;
my $ss = VRPipe::StepState->create(stepmember => 1, dataelement => 1, pipelinesetup => 1);
my @job_args;
my @sub_args;

for my $i (1 .. 10) {
    push(@job_args, { cmd => "job $i", dir => '/tmp' });
    push(@sub_args, { job => $i, stepstate => 1, requirements => 1, scheduler => 1 });
}
VRPipe::Job->bulk_create_or_update(@job_args);
VRPipe::Submission->bulk_create_or_update(@sub_args);

my @submissions = VRPipe::Submission->search({});
is_deeply [map { $_->id } @submissions], [1 .. 10], 'created 10 submissions to test with';

my $fm          = Parallel::ForkManager->new(4);
my $num_claimed = 0;
$fm->run_on_finish(
    sub {
        my ($pid, $claimed) = @_;
        $num_claimed++ if $claimed;
    }
);
for (1 .. 4) {
    $fm->start and next;
    
    my $claimed = $submissions[0]->_get_claim;
    
    $fm->finish($claimed);
}
$fm->wait_all_children;

is $num_claimed, 1, 'when we attempt to _get_claim the same submission simultaneously in 4 processes, only 1 succeeds';
$submissions[0]->reselect_values_from_db;
is $submissions[0]->_claim, 1, 'the submission object has stored the claim on itself';

$num_claimed = 0;
for (1 .. 4) {
    $fm->start and next;
    
    my $claimed = $submissions[0]->_get_claim;
    
    $fm->finish($claimed);
}
$fm->wait_all_children;

is $num_claimed, 0, 'when we attempt to _get_claim the same submission again, none succeed';

# do something more like vrpipe-handler does
$num_claimed = 0;
for (1 .. 4) {
    $fm->start and next;
    
    my $claimed_and_ran = $submissions[1]->claim_and_run(allowed_time => 5);
    
    $fm->finish($claimed_and_ran);
}
$fm->wait_all_children;

is $num_claimed, 1, 'when we attempt to claim_and_run the same submission simultaneously in 4 processes, exactly 1 succeeds';
$submissions[1]->reselect_values_from_db;
my $run_job = $submissions[1]->job;
is_deeply [$submissions[1]->_claim, defined $run_job->start_time], [1, 1], "claim_and_run resulted in the submission's _claim being true and the job's start_time being set";

# make sure that if a process dies in such a way that a submission is left with
# _claim true and the job dead, we can re-claim the submission and try to run
# the job again
$run_job->beat_heart;
sleep($run_job->survival_time + 1);

$num_claimed = 0;
for (1 .. 4) {
    $fm->start and next;
    
    my $claimed_and_ran = $submissions[1]->claim_and_run(allowed_time => 5);
    
    $fm->finish($claimed_and_ran);
}
$fm->wait_all_children;

is $num_claimed, 1, 'we were able to reclaim a submission when the _claim got stuck on and the job died';

$run_job->reselect_values_from_db;
$run_job->beat_heart;
sleep($run_job->survival_time + 1);
$run_job->beat_heart;

$num_claimed = 0;
for (1 .. 4) {
    $fm->start and next;
    
    my $claimed_and_ran = $submissions[1]->claim_and_run(allowed_time => 5);
    
    $fm->finish($claimed_and_ran);
}
$fm->wait_all_children;

is $num_claimed, 0, 'but we do not reclaim a submission when it is running ok';

exit;
