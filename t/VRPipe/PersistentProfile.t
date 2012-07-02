#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use TryCatch;

BEGIN {
    use Test::Most tests => 4;
    # only the author needs to run this test
    use VRPipeTest (required_env => 'VRPIPE_TEST_PIPELINES');
    use TestPipelines;
    
    use_ok('VRPipe::Persistent::Schema');
}

my ($l1, $t1);

my $j1 = VRPipe::Job->get(cmd => 'foo');
my $j2 = VRPipe::Job->get(cmd => 'foo');
is $j2->id, $j1->id, 'create followed by find worked';

start_clock(__LINE__);
for my $i (1..1000) {
    VRPipe::Job->get(cmd => qq[echo "job $i ]. 'n'x1000 .qq[";]);
}
lap(__LINE__);

start_clock(__LINE__);
for my $i (1..1000) {
    VRPipe::Job->get(cmd => qq[echo "job $i ]. 'n'x1000 .qq[";]);
}
lap(__LINE__);

start_clock(__LINE__);
for my $i (1..1000) {
    my $file = VRPipe::File->get(path => "/my/file/path/$i");
    $file->add_metadata({foo => 'bar',
                         baz => 'loman',
                         cat => 'dog',
                         fish => 'turnip'});
}
lap(__LINE__);

start_clock(__LINE__);
for my $i (1..1000) {
    my $file = VRPipe::File->get(path => "/my/file/path/$i");
    $file->add_metadata({foo => 'boo',
                         rat => 'king'.$i});
}
lap(__LINE__);

start_clock(__LINE__);
my $correct_meta = 0;
for my $i (1..1000) {
    my $file = VRPipe::File->get(path => "/my/file/path/$i");
    my $metadata = $file->metadata;
    $correct_meta++ if ($metadata->{rat} eq 'king'.$i && $metadata->{foo} eq 'boo' && $metadata->{baz} eq 'loman');
}
lap(__LINE__);
is $correct_meta, 1000, 'basic file metadata test worked';

start_clock(__LINE__);
for my $i (1..1000) {
    my $file = VRPipe::File->get(path => "/my/file/path/$i");
}
lap(__LINE__);

# do the critical bits that Manager currently does to submit jobs; first set
# up the objects
start_clock(__LINE__);
my $scheduler = VRPipe::Scheduler->get;
my @reqs = map { VRPipe::Requirements->get(memory => $_->[0], time => $_->[1]) } ([500, 1], [1000, 2], [2000, 3]);
my $pipeline = VRPipe::Pipeline->get(name => 'test_pipeline');
$pipeline->steps;
my $ds = VRPipe::DataSource->get(type => 'list', method => 'all', source => file(qw(t data datasource.onelist)));
my $ps = VRPipe::PipelineSetup->get(name => 'ps', datasource => $ds, output_root => dir(qw(tmp)), pipeline => $pipeline, options => {});
$ds->elements;
my $ss = VRPipe::StepState->get(stepmember => 1, dataelement => 1, pipelinesetup => 1);
my @job_args;
my @sub_args;
my $job_offset = 1001;
for my $i (1..10000) {
    push(@job_args, { cmd => "job $i", dir => '/tmp' });
    push(@sub_args, { job => $i + $job_offset, stepstate => 1, requirements => 1, scheduler => 1, '_done' => 1 });
}
for my $i (10001..11000) {
    push(@job_args, { cmd => "job $i", dir => '/tmp', $i > 10980 ? (running => 1) : () });
    push(@sub_args, { job => $i + $job_offset, stepstate => 1, requirements => $i % 2 ? 2 : 3, scheduler => 1, '_done' => 0, $i > 10980 ? ('_sid' => 20) : () });
}
for my $i (11001..12000) {
    push(@job_args, { cmd => "job $i", dir => '/tmp' });
    push(@sub_args, { job => $i + $job_offset, stepstate => 1, requirements => $i % 2 ? 3 : 2, scheduler => 1, '_done' => 0, $i > 11980 ? ('_failed' => 1) : () });
}
VRPipe::Job->bulk_create_or_update(@job_args);
VRPipe::Submission->bulk_create_or_update(@sub_args);
lap(__LINE__);

# now do something like Manager and Scheduler do
start_clock(__LINE__);
my $added = 0;
my $count = VRPipe::Job->search({ 'running' => 1 });
my ($last_sub_id) = VRPipe::Submission->get_column_values('id', {}, { order_by => { -desc => 'id' }, rows => 1 });
my $pager = VRPipe::Submission->search_paged({ '_done' => 0, '_failed' => 0, '_sid' => undef, 'id' => { '<=' => $last_sub_id } }, { order_by => 'requirements' }, 1000);
while (my $subs = $pager->next) {
    my %batches;
    foreach my $sub (@$subs) {
        #next if $sub->sid;
        my $job = $sub->job;
        next if $job->block_and_skip_if_ok;
        push(@{$batches{$sub->requirements->id}}, $sub);
    }
    
    while (my ($req_id, $batch_subs) = each %batches) {
        submit($batch_subs, $batch_subs->[0]->requirements);
    }
    
    if ($added == 0) {
        my $job = VRPipe::Job->get(cmd => "job extra", dir => '/tmp');
        VRPipe::Submission->get(job => $job, stepstate => 1, requirements => 3, scheduler => 1, '_done' => 0, '_failed' => 0);
        
        my ($submission) = VRPipe::Submission->search({ '_failed' => 1, requirements => 2 }, { rows => 1 });
        $submission->_failed(0);
        $submission->update;
        
        $added = 1;
    }
}
sub submit {
    my ($submissions, $requirements) = @_;
    
    my $aid = 0;
    my $parray = VRPipe::PersistentArray->get(members => $submissions);
    my $for = $parray;
    
    my $all_claimed = 1;
    foreach my $sub (@$submissions) {
        my $claimed = $sub->claim;
        unless ($claimed) {
            $all_claimed = 0;
            last;
        }
        $sub->_hid($for->id);
        $sub->_aid(++$aid);
    }
    
    my $sid = $requirements->id.'1234';
    if ($all_claimed) {
        foreach my $sub (@$submissions) {
            $sub->sid($sid);
        }
    }
}
is_deeply [scalar(VRPipe::Submission->search({ '_sid' => 21234 })), scalar(VRPipe::Submission->search({ '_sid' => 31234 }))], [981, 980], 'all submissions were submitted';
lap(__LINE__);

done_testing;
exit;

sub start_clock {
    $l1 = shift;
    $t1 = time();
}

sub lap {
    my $l2 = shift;
    my $t2 = time();
    note("Going from line $l1..$l2 took ", $t2 - $t1, " seconds\n");
    $t1 = time();
    $l1 = $l2 + 1;
}