#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 3;
    
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