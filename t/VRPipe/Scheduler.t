#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Sys::Hostname;

BEGIN {
    use Test::Most tests => 11;
    use VRPipeTest;
    
    use_ok('VRPipe::Scheduler');
}

ok my $scheduler = VRPipe::Scheduler->create, 'able to create a VRPipe::Scheduler based on configured default';

# local
ok $scheduler = VRPipe::Scheduler->create(type => 'local'), q[able to get the lsf scheduler using get(type => 'local')];
is $scheduler->type, 'local', 'the type really is local';

my $requirements = VRPipe::Requirements->create(memory => 1, time => 1);
is $scheduler->determine_queue($requirements), 'local', q[determine_queue() returned local's only queue];

# lsf
ok $scheduler = VRPipe::Scheduler->create(type => 'lsf'), q[able to get the lsf scheduler using get(type => 'lsf')];
is $scheduler->type, 'lsf', 'the type really is lsf';

SKIP: {
    my $host = hostname();
    skip "author-only lsf tests", 4 unless $host eq 'uk10k-1-1-01';
    
    is $scheduler->determine_queue($requirements), 'normal', 'determine_queue() gave normal queue for 10MB and 1hr';
    $requirements = VRPipe::Requirements->create(memory => 37000, time => 1);
    is $scheduler->determine_queue($requirements), 'test', 'determine_queue() gave test queue for 37GB and 1hr';
    $requirements = VRPipe::Requirements->create(memory => 1, time => 13);
    is $scheduler->determine_queue($requirements), 'long', 'determine_queue() gave long queue for 10MB and 13hr';
    $requirements = VRPipe::Requirements->create(memory => 1, time => 49);
    is $scheduler->determine_queue($requirements), 'basement', 'determine_queue() gave basement queue for 10MB and 49hr';
}

done_testing;
exit;
