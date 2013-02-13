#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Sys::Hostname;

BEGIN {
    use Test::Most tests => 15;
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
    is $scheduler->determine_queue($requirements), 'normal', 'determine_queue() gave test queue for 37GB and 1hr'; # used to be 'test' before our memory limits were removed from all queues
    $requirements = VRPipe::Requirements->create(memory => 1, time => 13);
    is $scheduler->determine_queue($requirements), 'long', 'determine_queue() gave long queue for 10MB and 13hr';
    $requirements = VRPipe::Requirements->create(memory => 1, time => 49);
    is $scheduler->determine_queue($requirements), 'basement', 'determine_queue() gave basement queue for 10MB and 49hr';
}

# ec2
SKIP: {
    eval "require VM::EC2;";
    skip "VM::EC2 is not installed", 4 if $@;
    
    ok $scheduler = VRPipe::Scheduler->create(type => 'ec2'), q[able to get the ec2 scheduler using get(type => 'ec2')];
    is $scheduler->type, 'ec2', 'the type really is ec2';
    
    $requirements = VRPipe::Requirements->create(memory => 100, time => 120);
    is $scheduler->determine_queue($requirements), 'm1.medium', 'determine_queue() gave m1.medium instance for 100MB and 2mins';
    
    my $scheduler_cmd_line = join(
        ' ',
        $scheduler->submit_command,
        $scheduler->submit_args(
            requirements => $requirements,
            stdo_file    => '/dev/null',
            stde_file    => '/dev/null',
            cmd          => 'the command to run'
        )
    );
    like $scheduler_cmd_line, qr/perl .+ -MVRPipe::Schedulers::ec2 -e "VRPipe::Schedulers::ec2->submit(@ARGV)" instance m1.medium memory 100 cmd 'the cmd to run'/, 'the expected scheduler cmd line could be constructed using submit_command() and submit_args()';
}

done_testing;
exit;
