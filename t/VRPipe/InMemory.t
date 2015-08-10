#!/usr/bin/env perl
use strict;
use warnings;
use Parallel::ForkManager;
use EV;
use AnyEvent;

BEGIN {
    use Test::Most tests => 143;
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

# test assert_life and is_alive
my $assert_pid = fork;
if (defined $assert_pid && $assert_pid == 0) {
    $im->assert_life('foo');
    sleep(2);
    exit(0);
}
sleep(1);
is $im->is_alive('foo'), 1, 'is_alive works while the pid that asserted life is alive';
waitpid($assert_pid, 0);
is $im->is_alive('foo'), 0, 'is_alive returns false once the pid that asserted life is dead';

ok $im->assert_life('bar'), 'assert_life() worked in parent pid';
is $im->is_alive('bar'),    1, 'is_alive() worked in the same process';
ok $im->assert_life('zar'), 'assert_life() worked in parent pid again, for a different key';
ok $im->assert_life('zar'), 'assert_life() succeeds with the same key on the same instance';
my $im2 = VRPipe::Persistent::InMemory->new;
ok !$im2->assert_life('zar'), 'assert_life() fails with the same key on a new instance';
ok $im2->assert_life('char'), 'assert_life() worked in a new instance, for a different key';
is $im2->is_alive('char'), 1, 'is_alive() worked on the new instance';
undef $im2;
is $im->is_alive('bar'),  1, 'after undefing the new instance, the orig instance assertions are still alive';
is $im->is_alive('char'), 0, 'but the new instance assertions are dead';

# test the log() method and its ability to email
warning_like { $im->log('log_test', email_admin => 1, force_when_testing => 1) } qr/^\d\d\d\d-\d\d-\d\d \d\d:\d\d:\d\d \| host \S+ \| pid \d+ \| log_test$/, 'log() generates a warning';
my @deliveries = Email::Sender::Simple->default_transport->deliveries;
my $email      = $deliveries[0]->{email};
my $email_body = $email->get_body;
$email_body =~ s/\s+$//;
is_deeply [scalar(@deliveries), $email->get_header("Subject"), $email_body], [1, 'VRPipe Server message', 'log_test'], 'log() also sent an email';
like $email->get_header('From'), qr/"VRPipe Server" <\S+\@\S+>/, 'from address looks good';

# check that we can log to redis
ok my $log_file = $VRPipe::Persistent::InMemory::log_file, 'class variable $log_file was set';
my $log_file_line = 0;
$im->log_stderr;
while ($redis->lpop('stderr')) {
    next; # clear out anything in there, ie. from the running vrpipe-server
}
latest_log_file_lines();
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

# test that we don't open zillions of Redis connections every time we call a
# Redis-backed method, but that we do create a new connection per process
my %redis_instances;
my $num_redis_instances = 0;
for my $i (1 .. 5) {
    my $redis = $im->_redis;
    unless (exists $redis_instances{"$redis"}) {
        $num_redis_instances++;
        $redis_instances{"$redis"} = 1;
    }
}
my $fm = Parallel::ForkManager->new(3);
$fm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my ($new_instances) = @$data_structure_reference;
        $num_redis_instances += $new_instances;
    }
);
for my $loop_num (1 .. 3) {
    $fm->start and next;
    my $new_instances = 0;
    my %own_instances;
    for my $i (1 .. 5) {
        my $redis = $im->_redis;
        unless (exists $redis_instances{"$redis"}) {
            unless (exists $own_instances{"$redis"}) {
                $new_instances++;
                $own_instances{"$redis"} = 1;
            }
        }
    }
    $fm->finish(0, [$new_instances]);
}
$fm->wait_all_children;
is $num_redis_instances, 4, 'new redis instances (connections) are only made when necessary';

# test locking
my $num_loops = 16;
$fm = Parallel::ForkManager->new($num_loops);
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

ok $im->lock('test_lock5', unlock_after => 4), 'lock attempt on new key works';
ok !$im->lock('test_lock5', unlock_after => 4), 'lock attempt after getting lock fails';
ok $im->locked('test_lock5'), 'locked() works immediately after locking';
sleep(2);
ok $im->locked('test_lock5'), 'locked() works 2s after locking';
sleep(3);
ok !$im->locked('test_lock5'), 'locked() returns false beyond unlock_after time';
ok $im->lock('test_lock5', unlock_after => 4), 'lock attempt on the expired key works';
$tested_positive = 0;
{
    $fm->start and next;
    my $locked = $im->locked('test_lock5');
    my $by_me  = $im->locked('test_lock5', by_me => 1);
    my $ok     = $locked && !$by_me;
    $fm->finish(0, [$ok]);
}
$fm->wait_all_children;
is $tested_positive, 1, 'locked(by_me) fails for another process while just locked() works';
ok $im->locked('test_lock5', by_me => 1), 'locked(by_me) works for the owner process';
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

$time = time();
ok $im->lock('test_lock6', unlock_after => 4), 'lock attempt on a new key worked';
ok $im->locked('test_lock6'), 'locked() still works immediately after locking';
$im->block_until_locked('test_lock6', check_every => 1);
$elapsed   = time() - $time;
$good_time = $elapsed <= 1;
ok $good_time, 'block_until_locked() did not make us wait on a key we just locked ourselves';
{
    $fm->start and next;
    $im->block_until_locked('test_lock6', check_every => 1);
    $fm->finish(0, []);
}
$fm->wait_all_children;
$elapsed = time() - $time;
$good_time = $elapsed >= 4 && $elapsed <= 6;
ok $good_time, 'block_until_locked() made us wait on a key locked by another process';

my $ev_timer = EV::timer 0, 0, sub {
    my $im = VRPipe::Persistent::InMemory->new;
    ok $im->lock('test_lock7', unlock_after => 2), 'lock attempt on a new key worked';
    my $cv = AnyEvent->condvar;
    my $sleep_timer = AnyEvent->timer(after => 3, cb => sub { $cv->send });
    $cv->recv;
    
    ok !$im->locked('test_lock7'), 'an lock became unlocked after waiting beyond initial expiry';
    
    $im->lock('test_lock7', unlock_after => 2);
    undef($im);
    $im = VRPipe::Persistent::InMemory->new;
    ok $im->locked('test_lock7'), 'a lock is still locked immediately after undefing the instance that created it';
    
    ok $im->unlock('test_lock7'), 'unlock() seemed to work';
    
    my $child_pid = fork();
    if (defined $child_pid && $child_pid == 0) {
        $im->lock('test_lock8', unlock_after => 2);
        exit(0);
    }
    waitpid($child_pid, 0);
    ok $im->locked('test_lock8'), 'a lock gained in a forked child with an expire time was kept';
    
    $cv = AnyEvent->condvar;
    $sleep_timer = AnyEvent->timer(after => 3, cb => sub { $cv->send });
    $cv->recv;
    
    ok !$im->locked('test_lock8'), 'and then lost again after the expire time';
    
    $child_pid = fork();
    if (defined $child_pid && $child_pid == 0) {
        $im->lock('test_lock9');
        exit(0);
    }
    waitpid($child_pid, 0);
    ok !$im->locked('test_lock9'), 'a lock gained in a forked child with no expire time was lost after the child exited';
    
    EV::unloop;
};
EV::run;

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

ok $im->note('test_note3', forget_after => 2, value => 'foo'), 'note(value) worked';
is $im->noted('test_note3'), 'foo', 'you can get the value back using noted()';
sleep(3);
ok !$im->noted('test_note3'), 'noted() returns false after the note expires, even when a value was supplied';

# briefly confirm that locking and notes work for Persistent class objects
my $file1 = VRPipe::File->create(path => '/foo');
my $file2 = VRPipe::File->create(path => '/bar');
my $req = VRPipe::Requirements->create(memory => 1, time => 1);
ok $file1->lock,   'was able to get a lock() on a VRPipe::Persistent instance with no args';
ok $file1->locked, 'locked() returns true for it';
$tested_positive = 0;
{
    $fm->start and next;
    my $locked = $file1->locked();
    my $by_me  = $file1->locked(by_me => 1);
    my $ok     = $locked && !$by_me;
    $fm->finish(0, [$ok]);
}
$fm->wait_all_children;
is $tested_positive, 1, 'locked(by_me) fails for another process while just locked() works';
ok $file1->locked(by_me => 1), 'locked(by_me) works for the owner process';
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

$time = time();
ok $req->lock(unlock_after => 4), 'lock attempt for a different Persistent class worked';
ok $req->locked, 'locked() returned true for it';
$req->block_until_unlocked(check_every => 1);
$elapsed = time() - $time;
ok $good_time, 'block_until_unlocked() made us wait until the lock expired';
ok !$req->locked, 'locked() returns false after waiting on the block';

$time = time();
ok $req->lock(unlock_after => 4), 'locked it again';
ok $req->locked, 'locked() returned true for it';
{
    $fm->start and next;
    $req->block_until_locked(check_every => 1, unlock_after => 2);
    $fm->finish(0, []);
}
$fm->wait_all_children;
$elapsed = time() - $time;
ok $good_time, 'block_until_locked() made us wait until the lock expired and was aquired by another process';
ok $req->locked, 'locked() returns true after waiting on the block';

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

# test sessions
ok my $session_key = $im->create_session({ foo => 'bar' }, idle_expiry => 2, max_life => 6), 'was able to use create_session()';
is_deeply $im->get_session($session_key), { foo => 'bar' }, 'get_session() worked prior to expiry';
ok $im->session_set($session_key, banana => 'rama'), 'session_set() could be called';
is_deeply $im->get_session($session_key), { foo => 'bar', banana => 'rama' }, 'session_set() really worked';
is $im->session_get($session_key, 'banana'), 'rama', 'session_get() works';
ok $im->session_del($session_key, 'foo'), 'session_del() could be called';
is_deeply $im->get_session($session_key), { banana => 'rama' }, 'session_del() really worked';
sleep(3);
is $im->get_session($session_key), undef, 'get_session() returned nothing after idle_expiry time';
ok my $session_key2 = $im->create_session({ cat => 'dog' }, idle_expiry => 2, max_life => 6), 'was able to create another session';
cmp_ok $session_key, 'ne', $session_key2, 'the second session key was not the same as the first';
my $got_session_ok = 0;

for (1 .. 7) {
    sleep(1);
    my $session = $im->get_session($session_key2);
    $got_session_ok++ if ($session && $session->{cat} eq 'dog');
}
is $got_session_ok, 6, 'calling get_session within idle_expiry time refreshed the expiry, but it eventually expired after max_life';

# test rate limiting
is $im->rate_limit('testlimit'), 1, 'rate_limit() returns true the first time';
is $im->rate_limit('testlimit'), 0, 'rate_limit() returns false the second time';
sleep(1);
is $im->rate_limit('testlimit'), 1, 'rate_limit() returns true after waiting a second';
is $im->rate_limit('testlimit2', punish_excess => 1), 1, 'rate_limit(punish_excess) returns true the first time';
is $im->rate_limit('testlimit2', punish_excess => 1), 0, 'rate_limit(punish_excess) returns false the second time';
sleep(1);
is $im->rate_limit('testlimit2', punish_excess => 1), 0, 'rate_limit(punish_excess) returns false the third time even after waiting a second';
sleep(3);
is $im->rate_limit('testlimit2', punish_excess => 1), 1, 'rate_limit(punish_excess) returns true after waiting 3 seconds';
is $im->rate_limit('testlimit3', per_second    => 2), 1, 'rate_limit(per_second) returns true the first time';
is $im->rate_limit('testlimit3', per_second    => 2), 1, 'rate_limit(per_second) returns true the second time';
is $im->rate_limit('testlimit3', per_second    => 2), 0, 'rate_limit(per_second) returns false the second time, since we only allowed 2 per second';
sleep(1);
is $im->rate_limit('testlimit3', per_second => 2), 1, 'rate_limit(per_second) returns true the fourth time after waiting a second';

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
