#!/usr/bin/env perl
use strict;
use warnings;
use Parallel::ForkManager;

BEGIN {
    use Test::Most tests => 98;
    use VRPipeTest;
    $ENV{EMAIL_SENDER_TRANSPORT} = 'Test';
    use_ok('VRPipe::Persistent::InMemory');
}

ok my $im = VRPipe::Persistent::InMemory->new, 'able to create an InMemory object';
isa_ok($im, 'VRPipe::Persistent::InMemory');

# check that we can get a redis object
ok my $redis = $im->_redis, 'able to call _redis()';
isa_ok($redis, 'Redis');
ok $redis->ping,      'Redis->ping worked';
ok $im->datastore_ok, 'datastore_ok() worked';

# test the log() method and its ability to email
warning_like { $im->log('log_test', email_admin => 1, force_when_testing => 1) } qr/^\d\d\d\d-\d\d-\d\d \d\d:\d\d:\d\d \| pid \d+ \| log_test$/, 'log() generates a warning';
my @deliveries = Email::Sender::Simple->default_transport->deliveries;
my $email      = $deliveries[0]->{email};
my $email_body = $email->get_body;
$email_body =~ s/\s+$//;
is_deeply [scalar(@deliveries), $email->get_header("Subject"), $email_body], [1, 'VRPipe Server message', 'log_test'], 'log() also sent an email';
like $email->get_header('From'), qr/"VRPipe Server" <\S+\@\S+>/, 'from address looks good';

# check that we can log to redis
ok my $log_file = $im->_log_file, '_log_file returning something';
my $log_file_line = 0;
latest_log_file_lines();
$im->log_stderr;
warn "This first warning should not appear on STDERR\n";
warn "This second warning should not appear on STDERR\n";
is $redis->lpop('stderr'), "This first warning should not appear on STDERR\n", 'the first warning was in redis';
warn "This third warning should not appear on STDERR\n";
$im->write_stderr;
is_deeply [latest_log_file_lines()], ["This second warning should not appear on STDERR\n", "This third warning should not appear on STDERR\n"], 'write_stderr wrote remaining warnings to the _log_file';

# test that log() works with redis, and also test debug()
$im->log('log to file');
$im->debug('debug without verbose');
is_deeply [latest_log_file_messages()], [], 'log() did not write directly to the log file';
$im->write_stderr;
is_deeply [latest_log_file_messages()], ['log to file'], 'log() wrote to the log file via redis, and debug() did nothing';
$im->verbose(1);
$im->debug('debug with verbose');
$im->write_stderr;
is_deeply [latest_log_file_messages()], ['debug with verbose'], 'debug() wrote to the log file via redis when verbose';
$im->verbose(0);

untie *STDERR;

# test locking
my $num_loops = 16;
my $fm        = Parallel::ForkManager->new($num_loops);
my $tested_positive;
$fm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my ($test_result) = @$data_structure_reference;
        $tested_positive++ if $test_result;
    }
);

$tested_positive = 0;
for my $loop_num (1 .. $num_loops) {
    $fm->start and next;
    my $got_lock = $im->lock('test_lock1', unlock_after => 3);
    $fm->finish(0, [$got_lock]);
}
$fm->wait_all_children;
is $tested_positive, 1, 'lock() on the same key 16 times in parallel only gave one process the lock';

$tested_positive = 0;
for my $loop_num (1 .. $num_loops) {
    $fm->start and next;
    my $got_lock = $im->lock('test_lock1', unlock_after => 3, non_exclusive => 1);
    $fm->finish(0, [$got_lock]);
}
$fm->wait_all_children;
is $tested_positive, 0, 'lock() on the same key 16 times in parallel with non_exclusive all failed since we had an exclusive lock';

$tested_positive = 0;
for my $loop_num (1 .. $num_loops) {
    $fm->start and next;
    my $got_lock = $im->lock('test_lock2', unlock_after => 2, non_exclusive => 1);
    $fm->finish(0, [$got_lock]);
}
$fm->wait_all_children;
is $tested_positive, 16, 'lock() on a new key 16 times in parallel with non_exclusive gave lock to all 16';

$tested_positive = 0;
for my $loop_num (1 .. $num_loops) {
    sleep(1) if $loop_num == 2;
    $fm->start and next;
    $im->lock('test_lock3', unlock_after => 5);
    sleep(3) if $loop_num == 1;
    my $unlocked = $im->unlock('test_lock3');
    $fm->finish(0, [$unlocked]);
}
$fm->wait_all_children;
is $tested_positive, 1, 'unlock() on the same key 16 times in parallel only worked on the process that got the lock';

$tested_positive = 0;
for my $loop_num (1 .. $num_loops) {
    sleep(1) if $loop_num == 2;
    $fm->start and next;
    $im->lock('test_lock4', unlock_after => 5);
    sleep(3) if $loop_num == 1;
    my $refreshed = $im->refresh_lock('test_lock4', unlock_after => 1);
    $fm->finish(0, [$refreshed]);
}
$fm->wait_all_children;
is $tested_positive, 1, 'refresh_lock() on the same key 16 times in parallel only worked on the process that got the lock';

ok $im->lock('test_lock5', unlock_after => 4), 'lock attempt on new key works';
ok !$im->lock('test_lock5', unlock_after => 4), 'lock attempt after getting lock fails';
ok $im->locked('test_lock5'), 'locked() works immediately after locking';
sleep(2);
ok $im->locked('test_lock5'), 'locked() works 2s after locking';
sleep(3);
ok !$im->locked('test_lock5'), 'locked() returns false beyond unlock_after time';
ok !$im->refresh_lock('test_lock5', unlock_after => 4), 'refresh_lock() fails after the lock expired';
ok $im->lock('test_lock5', unlock_after => 4), 'lock attempt on the expired key works';
sleep(2);
ok $im->locked('test_lock5'), 'locked() works again 2s after locking';
ok $im->refresh_lock('test_lock5', unlock_after => 4), 'refresh_lock() works while locked';
sleep(3);
ok $im->locked('test_lock5'), 'locked() works after initial expire time thanks to refresh updating the timeout';
sleep(2);
ok !$im->locked('test_lock5'), 'locked() returns false beyond the refresh time';
ok $im->lock('test_lock5', unlock_after => 4), 'lock attempt on the expired key works again';
ok $im->locked('test_lock5'), 'locked() still works immediately after locking';
ok $im->unlock('test_lock5'), 'unlock() seems to have worked';
ok !$im->locked('test_lock5'), 'locked() returns false after unlocking';

my $time = time();
ok $im->lock('test_lock5', unlock_after => 4), 'lock attempt on the expired key works again';
ok $im->locked('test_lock5'), 'locked() still works immediately after locking';
$im->block_until_unlocked('test_lock5', check_every => 1);
my $elapsed = time() - $time;
my $good_time = $elapsed == 4 || $elapsed == 5;
ok $good_time, 'block_until_unlocked() made us wait until the lock expired';
ok !$im->locked('test_lock5'), 'locked() returns false after waiting on the block';

# test note taking
$tested_positive = 0;
for my $loop_num (1 .. $num_loops) {
    $fm->start and next;
    my $set_note = $im->note('test_note1', forget_after => 2);
    $fm->finish(0, [$set_note]);
}
$fm->wait_all_children;
is $tested_positive, 16, 'note() on the same key 16 times in parallel all worked';

ok $im->note('test_note2', forget_after => 2), 'note() on a new key worked';
ok $im->noted('test_note2'), 'noted() works immediately after setting the note';
sleep(3);
ok !$im->noted('test_note2'), 'noted() returns false after the note expires';
ok $im->note('test_note2', forget_after => 2), 'note() on the same key worked again';
sleep(1);
ok $im->noted('test_note2'), 'noted() works 1s after setting the note';
ok $im->note('test_note2', forget_after => 4), 'note() on the same key worked again while the key was still active';
sleep(3);
ok $im->noted('test_note2'), 'noted() works after the initial expire time of the note';
sleep(2);
ok !$im->noted('test_note2'), 'noted() returns false after the note expires again';
ok $im->note('test_note2', forget_after => 2), 'note() on the same key worked again';
ok $im->noted('test_note2'),       'noted() works immediately after setting the note';
ok $im->forget_note('test_note2'), 'forget_note() seemed to work';
ok !$im->noted('test_note2'), 'noted() returns false after using forget_note()';

# briefly confirm that locking and notes work for Persistent class objects
my $file1 = VRPipe::File->create(path => '/foo');
my $file2 = VRPipe::File->create(path => '/bar');
my $req = VRPipe::Requirements->create(memory => 1, time => 1);
ok $file1->lock,   'was able to get a lock() on a VRPipe::Persistent instance with no args';
ok $file1->locked, 'locked() returns true for it';
ok !$file2->locked, 'locked() returns false for an unlocked Persistent instance with a different id';
ok !$req->locked,   'locked() returns false for an unlocked Persistent instance with the same id but different class';
ok $file1->unlock, 'unlock() worked on the locked Persistent';
ok !$file1->locked, 'locked() returns false after unlocking';
ok $file2->lock(unlock_after => 2), 'was able to lock the other Persistent instance, passing in unlock_after arg';
ok $file2->locked, 'immediately after locking it it was locked()';
sleep(3);
ok !$file2->locked, 'it became unlocked after the unlock_after time';
ok $file2->lock(unlock_after => 2), 'was able to lock the other Persistent instance again';
sleep(1);
ok $file2->locked, '1s after locking it it was locked()';
ok $file2->refresh_lock(unlock_after => 3), 'refresh_lock() seemed to work';
sleep(2);
ok $file2->locked, 'still locked after initial lock time';
sleep(2);
ok !$file2->locked, 'it became unlocked after the refreshed time';

$time = time();
ok $req->lock(unlock_after => 4), 'lock attempt for a different Persistent class worked';
ok $req->locked, 'locked() returned true for it';
$req->block_until_unlocked(check_every => 1);
$elapsed = time() - $time;
ok $good_time, 'block_until_unlocked() made us wait until the lock expired';
ok !$req->locked, 'locked() returns false after waiting on the block';

ok $file1->note('test_note'),  'was able to set a note() on a VRPipe::Persistent instance';
ok $file1->noted('test_note'), 'noted() returns true for it';
ok !$file2->noted('test_note'), 'noted() returns false for an un-noted Persistent instance with a different id';
ok !$req->noted('test_note'),   'noted() returns false for an un-noted Persistent instance with the same id but different class';
ok $file1->forget_note('test_note'), 'forget_note() seemed to work';
ok !$file1->noted('test_note'), 'noted() returns false after forgetting';

ok $file1->lock(unlock_after => 2), 'was able to get a lock() on a VRPipe::Persistent instance again';
ok $file1->note('test_note2', forget_after => 2), 'was able to set a note() on a locked instance';
ok $file1->locked, 'locked() returns true';
ok $file1->noted('test_note2'), 'noted() returns true';

# test queuing
$redis->del('lock.test_queue');
$redis->del('queue.test_queue');
$redis->del('queue.test_queue2');
ok $im->lock('test_queue'), 'was able to get a lock with name test_queue';
is_deeply [$im->queue('test_queue')], [], 'queue() returns empty array for queue test_queue';
ok $im->enqueue('test_queue', 'foo'), 'enqueue() seemed to work';
ok $im->unlock('test_queue'), 'unlock test_queue worked; if subsequent test works we know locks and queues are independent';
is_deeply [$im->queue('test_queue')], ['foo'], 'queue() returns the value we enqueued';
ok $im->enqueue('test_queue',  'bar'),   'enqueue() worked a second time for the same queue';
ok $im->enqueue('test_queue2', 'lemon'), 'enqueue() worked for a new queue';
is_deeply [sort $im->queue('test_queue')], ['bar', 'foo'], 'queue() returns the values we enqueued';
is_deeply [$im->queue('test_queue2')], ['lemon'], 'queue() returns the values we enqueued';
ok $im->dequeue('test_queue', ['foo']), 'dequeue() with a value seemed to work';
is_deeply [$im->queue('test_queue')], ['bar'], 'queue() returns the values we enqueued minus the one we dequeued';
ok $im->drop_queue('test_queue'), 'drop_queue() seemed to work';
is_deeply [$im->queue('test_queue')], [], 'queue() returns empty array after dropping';
ok $im->enqueue('test_queue2', 'sand'), 'enqueue() worked again';
ok my $dequeued = $im->dequeue('test_queue2'), 'dequeue() with no value returned something';
my %allowed = (lemon => 1, sand => 1);
ok exists $allowed{$dequeued}, 'the dequeued value we got was one of those that were enqueued';
delete $allowed{$dequeued};
is_deeply [$im->queue('test_queue2')], [sort keys %allowed], 'queue() returns the remaining value we had enqueued';

done_testing;
exit;

sub latest_log_file_lines {
    open(my $fh, $log_file) || return;
    my @std_err;
    my $this_line = 0;
    while (<$fh>) {
        $this_line++;
        if ($this_line > $log_file_line) {
            push(@std_err, $_);
        }
    }
    $log_file_line = $this_line;
    return @std_err;
}

sub latest_log_file_messages {
    my @msgs;
    foreach my $line (latest_log_file_lines()) {
        my ($msg) = $line =~ /\| pid \d+ \| (.+)/;
        push(@msgs, $msg);
    }
    return @msgs;
}