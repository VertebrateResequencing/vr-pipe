#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 15;
    use VRPipeTest;
}

my $subject = VRPipe::MessageTracker->digest('overall setup state for setup 1');
is $subject, 'dc2e77c241e9a263a56eb1111817f038', 'digest() returned the correct value';

my $mt = VRPipe::MessageTracker->create(subject => $subject);
is $mt->id, 1, 'we were able to create a new MessageTracker object in the db';
is $mt->subject, $subject, 'it has the subject we set';
is $mt->message, '', 'it starts off with an empty message field';

is $mt->already_sent('complete with 100 elements'), 0, 'already_sent with a new message returns 0';
is $mt->message, '601dcdf8d4318546d86e745c8531c586', 'already_sent actually changed the stored message';
is $mt->already_sent('complete with 100 elements'), 1, 'already_sent with the same message returns 1';
is $mt->already_sent('complete with 101 elements'), 0, 'already_sent with a new message returns 0 again';
is $mt->already_sent('complete with 101 elements'), 1, 'already_sent with the same message returns 1 again';

my $mt2 = VRPipe::MessageTracker->create(subject => 'overall setup state for setup 2');
is $mt2->id, 2, 'we were able to create a second MessageTracker object in the db, passing in a raw subject to create()';
is $mt2->subject, '9438f03d0306250d65142c789b7787cc', 'the subject stored in the db is the digested form though';

my $mt3 = VRPipe::MessageTracker->get(subject => 'overall setup state for setup 1');
is_deeply [$mt3->id, $mt3->subject], [1, $subject], 'we can get() with the raw subject string as well';

is $mt3->message, '65eba99f537b1391741560626b07e673', 'message was as expected before update_message()';
$mt3->update_message('foobar');
is $mt3->message, '3858f62230ac3c915f300c664312c63f', 'update_message() changed message in db';
is $mt3->already_sent('foobar'), 1, 'passing the updated message to already_sent makes it return 1';

exit;
