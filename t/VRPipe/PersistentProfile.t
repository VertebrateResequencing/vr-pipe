#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use TryCatch;
use Time::HiRes qw(gettimeofday tv_interval);

BEGIN {
    use Test::Most tests => 3;
    # only the author needs to run this test
    use VRPipeTest (required_env => 'VRPIPE_TEST_PIPELINES');
    use TestPipelines;
    
    use_ok('VRPipe::Persistent::Schema');
}

my %times;
my %elapsed;

my $j1 = VRPipe::Job->create(cmd => 'foo');
my $j2 = VRPipe::Job->get(cmd => 'foo');
is $j2->id, $j1->id, 'create followed by get worked';

my $l = start_clock(__LINE__);
for my $i (1 .. 1000) {
    VRPipe::Job->create(cmd => qq[echo "job $i ] . 'n' x 1000 . qq[";]);
}
elapsed($l, __LINE__);

$l = start_clock(__LINE__);
for my $i (1 .. 1000) {
    VRPipe::Job->get(cmd => qq[echo "job $i ] . 'n' x 1000 . qq[";]);
}
elapsed($l, __LINE__);

$l = start_clock(__LINE__);
for my $i (1 .. 1000) {
    my $file = VRPipe::File->create(path => "/my/file/path/$i");
    $file->add_metadata({
            foo  => 'bar',
            baz  => 'loman',
            cat  => 'dog',
            fish => 'turnip'
        }
    );
}
elapsed($l, __LINE__);

$l = start_clock(__LINE__);
for my $i (1 .. 1000) {
    my $file = VRPipe::File->get(path => "/my/file/path/$i");
    $file->add_metadata({
            foo => 'boo',
            rat => 'king' . $i
        }
    );
}
elapsed($l, __LINE__);

$l = start_clock(__LINE__);
my $correct_meta = 0;
for my $i (1 .. 1000) {
    my $file = VRPipe::File->get(path => "/my/file/path/$i");
    my $metadata = $file->metadata;
    $correct_meta++ if ($metadata->{rat} eq 'king' . $i && $metadata->{foo} eq 'boo' && $metadata->{baz} eq 'loman');
}
elapsed($l, __LINE__);
is $correct_meta, 1000, 'basic file metadata test worked';

$l = start_clock(__LINE__);
for my $i (1 .. 1000) {
    my $file = VRPipe::File->get(path => "/my/file/path/$i");
}
elapsed($l, __LINE__);

# do the critical bits that Manager currently does to submit jobs; first set
# up the objects
$l = start_clock(__LINE__);
my $scheduler = VRPipe::Scheduler->create;
my @reqs      = map { VRPipe::Requirements->create(memory => $_->[0], time => $_->[1]) } ([500, 1], [1000, 2], [2000, 3]);
my $pipeline  = VRPipe::Pipeline->create(name => 'test_pipeline');
my $ds        = VRPipe::DataSource->create(type => 'list', method => 'all', source => file(qw(t data datasource.onelist))->absolute);
my $ps        = VRPipe::PipelineSetup->create(name => 'ps', datasource => $ds, output_root => dir(qw(tmp)), pipeline => $pipeline, options => {});
$ps->active(0);
$ps->update;
$ds->elements;
my $ss = VRPipe::StepState->create(stepmember => 1, dataelement => 1, pipelinesetup => 1);
my @job_args;
my @sub_args;
my $job_offset = 1001;

for my $i (1 .. 10000) {
    push(@job_args, { cmd => "job $i", dir => '/tmp' });
    push(@sub_args, { job => $i + $job_offset, stepstate => 1, requirements => 1, scheduler => 1, '_done' => 1 });
}
for my $i (10001 .. 11000) {
    push(@job_args, { cmd => "job $i", dir => '/tmp', $i > 10980 ? (start_time => DateTime->now, heartbeat => DateTime->now) : () });
    push(@sub_args, { job => $i + $job_offset, stepstate => 1, requirements => $i % 2 ? 2 : 3, scheduler => 1, '_done' => 0, $i > 10980 ? ('_claim' => 1) : () });
}
for my $i (11001 .. 12000) {
    push(@job_args, { cmd => "job $i", dir => '/tmp' });
    push(@sub_args, { job => $i + $job_offset, stepstate => 1, requirements => $i % 2 ? 3 : 2, scheduler => 1, '_done' => 0, $i > 11980 ? ('_failed' => 1) : () });
}
VRPipe::Job->bulk_create_or_update(@job_args);
VRPipe::Submission->bulk_create_or_update(@sub_args);
elapsed($l, __LINE__);

# now do something like Manager and Scheduler used to do *** these tests need
# updating to do something like vrpipe-handler does instead
#*** and they're disabled now because they're too slow, due to 1-second hack
# that ensures claim_and_run only succeeds for 1 process at a time; in the real
# world this is ok since claim_and_run is called in the farm, not in a single-
# process loop
# $l = start_clock(__LINE__);
# my $added = 0;
# my $count = VRPipe::Job->search({ 'heartbeat' => { '!=' => undef } });
# my ($last_sub_id) = VRPipe::Submission->get_column_values('id', {}, { order_by => { -desc => 'id' }, rows => 1 });
# my $pager = VRPipe::Submission->search_paged({ '_done' => 0, '_failed' => 0, '_claim' => 0, 'me.id' => { '<=' => $last_sub_id } }, { order_by => 'requirements', prefetch => [qw(job requirements)] }, 1000);
# while (my $subs = $pager->next) {
#     my $sub_loop = start_clock(__LINE__);
#     my %batches;
#     foreach my $sub (@$subs) {
#         my $job = $sub->job;
#         next if $job->block_and_skip_if_ok;
#         push(@{ $batches{ $sub->requirements->id } }, $sub);
#     }
#     elapsed($sub_loop, __LINE__, 1);

#     my $batch_loop = start_clock(__LINE__);
#     while (my ($req_id, $batch_subs) = each %batches) {
#         my $submit_call = start_clock(__LINE__);
#         submit($batch_subs, $batch_subs->[0]->requirements);
#         elapsed($submit_call, __LINE__, 1);
#     }
#     elapsed($batch_loop, __LINE__, 1);

#     if ($added == 0) {
#         my $job = VRPipe::Job->create(cmd => "job extra", dir => '/tmp');
#         VRPipe::Submission->create(job => $job, stepstate => 1, requirements => 3, scheduler => 1, '_done' => 0, '_failed' => 0);

#         my ($submission) = VRPipe::Submission->search({ '_failed' => 1, requirements => 2 }, { rows => 1 });
#         $submission->_failed(0);
#         $submission->update;

#         $added = 1;
#     }
# # }

# sub submit {
#     my ($submissions, $requirements) = @_;

#     my $create_call = start_clock(__LINE__);
#     my $parray = VRPipe::PersistentArray->create(members => $submissions);
#     elapsed($create_call, __LINE__, 1);
#     my $for = $parray;

#     my $aid             = 0;
#     my $all_claimed     = 1;
#     my $second_sub_loop = start_clock(__LINE__);
#     foreach my $sub (@$submissions) {
#         my ($claimed) = $sub->claim_and_run;
#         unless ($claimed == 1) {
#             $all_claimed = 0;
#             last;
#         }
#     }
#     elapsed($second_sub_loop, __LINE__, 1);
# }
# is scalar(VRPipe::Submission->search({ '_claim' => 1 })), 1981, 'all submissions were claimed';
# elapsed($l, __LINE__);

report();

done_testing;
exit;

sub start_clock {
    my $l1 = shift;
    $times{$l1} = [gettimeofday];
    return $l1;
}

sub elapsed {
    my ($l1, $l2, $silent) = @_;
    my $e = sprintf("%0.2f", tv_interval($times{$l1}));
    my $id = "$l1..$l2";
    note("Going from line $id took $e seconds\n") unless $silent;
    push(@{ $elapsed{$id} }, $e);
}

sub report {
    my %cums;
    while (my ($id, $times) = each %elapsed) {
        my $cum = 0;
        foreach my $t (@$times) {
            $cum += $t;
        }
        $cums{$id} = $cum;
    }
    
    note("\nMost time consuming sections:\n");
    foreach my $id (sort { $cums{$b} <=> $cums{$a} || $a cmp $b } keys %cums) {
        my $cum = $cums{$id};
        
        my $times = $elapsed{$id};
        my $total = 0;
        my $count = 0;
        foreach my $t (@$times) {
            $count++;
            $total += $t;
        }
        my $avg = sprintf("%0.2f", $total / $count);
        
        my $note = "\t$id: $cum seconds";
        if ($count > 1) {
            $note .= " ($avg avg over $count loops)\n";
        }
        else {
            $note .= "\n";
        }
        note($note);
    }
}
