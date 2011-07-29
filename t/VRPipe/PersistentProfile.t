#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 2;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPersistentReal;
}

my $j1 = VRPipe::Job->get(cmd => 'foo');
my $j2 = VRPipe::Job->get(cmd => 'foo');
is $j2->id, $j1->id, 'create followed by find worked';

my ($l1, $t1);

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

done_testing;
exit;

sub start_clock {
    $l1 = shift;
    $t1 = time();
}

sub lap {
    my $l2 = shift;
    my $t2 = time();
    warn "Going from line $l1..$l2 took ", $t2 - $t1, " seconds\n";
    $t1 = time();
    $l1 = $l2 + 1;
}