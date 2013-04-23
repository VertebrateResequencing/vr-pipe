#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class qw(file dir);
use Parallel::ForkManager;

BEGIN {
    use Test::Most tests => 8;
    use VRPipeTest;
    use TestPipelines;
}

# make a stepstate for some basic testing
my $test_pipeline = VRPipe::Pipeline->create(name => 'test_pipeline');
my $ds = VRPipe::DataSource->create(type => 'list', method => 'all', source => file(qw(t data datasource.onelist))->absolute);
my $ps = VRPipe::PipelineSetup->create(name => 'ps', datasource => $ds, output_root => dir(qw(tmp)), pipeline => $test_pipeline, options => {});
$ps->active(0); # stop the server trying to do something with this ps
$ps->update;
$ds->elements;
my $ss = VRPipe::StepState->create(stepmember => 1, dataelement => 1, pipelinesetup => 1);
my $des = VRPipe::DataElementState->get(id => 1);

# first just check that the basics of how we expect transactions and row locking
# to work
my $fm                 = Parallel::ForkManager->new(3);
my $did_complete       = 0;
my $already_complete   = 0;
my $confirmed_complete = 0;
$fm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my ($dc, $ac, $cc) = @$data_structure_reference;
        $did_complete++       if $dc;
        $already_complete++   if $ac;
        $confirmed_complete++ if $cc;
    }
);
for (1 .. 3) {
    $did_complete       = 0;
    $already_complete   = 0;
    $confirmed_complete = 0;
    $ss->reselect_values_from_db;
    $ss->complete(0);
    $ss->update;
    $des->reselect_values_from_db;
    $des->completed_steps(0);
    $des->update;
    
    foreach my $i (1 .. 3) {
        $fm->start and next;
        
        #sleep($i);
        
        my ($dc, $ac, $cc);
        #warn "$$ loop start\n";
        my $transaction = sub {
            #warn "$$ transaction start\n";
            #my ($state) = VRPipe::StepState->search({ id => $ss->id }, { for => 'update' });
            $ss->lock_row($ss);
            #$state->reselect_values_from_db;
            
            if ($ss->complete) {
                #warn "$$ already complete\n";
                $ac = 1;
            }
            else {
                sleep(5);
                #warn "$$ completing\n";
                if (nested_transaction()) {
                    $des->reselect_values_from_db;
                    if ($des->completed_steps == 1) {
                        $ss->complete(1);
                        $ss->update;
                        $dc = 1;
                        
                        if (nested_transaction_two()) {
                            $cc = 1;
                        }
                    }
                }
            }
            #warn "$$ transaction end\n";
        };
        $ss->do_transaction($transaction, 'failed');
        #warn "$$ loop end\n";
        $fm->finish(0, [$dc, $ac, $cc]);
    }
    $fm->wait_all_children;
    is_deeply [$did_complete, $already_complete, $confirmed_complete], [1, 2, 1], 'basics of transactions and row locking work';
}

sub nested_transaction {
    my $at_zero     = 0;
    my $transaction = sub {
        my $des = VRPipe::DataElementState->get(id => 1);
        $ss->lock_row($des);
        if ($des->completed_steps == 0) {
            $at_zero = 1;
            $des->completed_steps(1);
            $des->update;
        }
    };
    $ss->do_transaction($transaction, 'failed nested');
    return $at_zero;
}

sub nested_transaction_two {
    my $found       = 0;
    my $transaction = sub {
        ($found) = VRPipe::StepState->get_column_values('id', { complete => 1 });
    };
    $ss->do_transaction($transaction, 'failed nested 2');
    return $found;
}

# now do a more real-world test that used to reveal an issue that broke the
# above locking assumptions and allowed multiple processes to work on the same
# stepstate at the same time, causing random havoc
my ($output_dir, $pipeline, $step) = create_single_step_pipeline('fake_fastq_metadata', 'fastq_files');
my $si_datasource = VRPipe::DataSource->create(
    type    => 'sequence_index',
    method  => 'lane_fastqs',
    source  => file(qw(t data datasource.sequence_index))->absolute,
    options => { local_root_dir => dir(".")->absolute->stringify }
);
my $setup = VRPipe::PipelineSetup->create(
    name        => 'fm_setup',
    datasource  => $si_datasource,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);

# make 2 more setups in different dirs; whilst we don't have any tests
# specific to these setups, creating them is sufficient to reveal problems
my $output_dir2 = get_output_dir('fm_setup2');
VRPipe::PipelineSetup->create(
    name        => 'fm_setup2',
    datasource  => $si_datasource,
    output_root => $output_dir2,
    pipeline    => $pipeline,
    options     => {}
);
my $output_dir3 = get_output_dir('fm_setup3');
VRPipe::PipelineSetup->create(
    name        => 'fm_setup3',
    datasource  => $si_datasource,
    output_root => $output_dir3,
    pipeline    => $pipeline,
    options     => {}
);

my @ofiles;
foreach my $basename (qw(2822_6.fastq 2822_6_1.fastq 2822_6_2.fastq 2822_7_1.fastq 2822_7_2.fastq 2823_4_1.fastq 2823_4_2.fastq 8324_8_1.fastq 8324_8_2.fastq)) {
    push(@ofiles, file('t', 'data', $basename)->absolute);
}

ok handle_pipeline(@ofiles), 'fastq_metadata pipeline ran ok and all input files still exist';
my $fastq = VRPipe::File->create(path => file(qw(t data 2822_6.fastq))->absolute);
is $fastq->metadata->{avg_read_length}, 10, 'it actually added some metadata';

# the main symptom of problems is that start_over happens
my @logs = VRPipe::PipelineSetupLog->search({ message => { like => "%StepState->start_over was called%" } });
is scalar(@logs), 0, 'no start_over calls were made';

# start_over, and other things going wrong, happen due to Jobs running more than
# once at a time, or triggers on the same stepstate happening more than once
# at a time. we don't mind the same pid repeating these things, since that
# suggests just a retry of a loop (serial repetition) instead of something
# actually bad happening (parallel repetition)
my @job_ids = VRPipe::Job->get_column_values('id', {});
my $ran_once_count = 0;
foreach my $jid (@job_ids) {
    my @logs = VRPipe::PipelineSetupLog->search({ job_id => $jid, message => { like => "%cmd-running child exited with code 0%" } });
    my %by_pid = map { $_->pid => 1 } @logs;
    my $count = keys %by_pid;
    $ran_once_count++ if $count == 1;
    if ($count != 1) {
        warn "job_id $jid is bad\n";
    }
}
is $ran_once_count, scalar(@job_ids), 'all jobs only ran once';

my @ss_ids = VRPipe::StepState->get_column_values('id', { pipelinesetup => { '!=' => 1 } });
my $triggered_once_count = 0;
foreach my $ss_id (@ss_ids) {
    my @logs = VRPipe::PipelineSetupLog->search({ ss_id => $ss_id, message => { like => "%so will trigger the next Step%" } });
    my %by_pid = map { $_->pid => 1 } @logs;
    my $count = keys %by_pid;
    $triggered_once_count++ if $count <= 1; # this should be == 1, but db issues mean we don't always create the Log rows we requested
    if ($count > 1) {                       # this should be != 1
        warn "ss_id $ss_id is bad\n";
    }
}
is $triggered_once_count, scalar(@ss_ids), 'all stepstates only triggered the next step once';

done_testing;
exit;
