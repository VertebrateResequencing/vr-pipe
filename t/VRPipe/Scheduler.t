#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Sys::Hostname;

BEGIN {
    use Test::Most tests => 20;
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
    skip "author-only lsf tests", 5 unless $host eq 'vr-2-2-02';
    
    is $scheduler->determine_queue($requirements), 'normal', 'determine_queue() gave normal queue for 10MB and 1hr';
    $requirements = VRPipe::Requirements->create(memory => 1, time => 300);
    is $scheduler->determine_queue($requirements), 'normal', 'determine_queue() gave normal queue for 10MB and 5mins';
    $requirements = VRPipe::Requirements->create(memory => 37000, time => 1);
    is $scheduler->determine_queue($requirements), 'normal', 'determine_queue() gave test queue for 37GB and 1hr'; # used to be 'test' before our memory limits were removed from all queues
    $requirements = VRPipe::Requirements->create(memory => 1, time => 13);
    is $scheduler->determine_queue($requirements), 'long', 'determine_queue() gave long queue for 10MB and 13hr';
    $requirements = VRPipe::Requirements->create(memory => 1, time => 262800);
    is $scheduler->determine_queue($requirements), 'basement', 'determine_queue() gave basement queue for 10MB and 73hr';
}

# ec2
SKIP: {
    eval "require VM::EC2;";
    skip "VM::EC2 is not installed", 5 if $@;
    # we might have VM::EC2 installed, but might not be using the ec2
    # scheduler, in which case it isn't configured and the module won't work
    my $vrp_config = VRPipe::Config->new();
    my $access_key = $vrp_config->ec2_access_key;
    skip "ec2 scheduler is not configured", 5 unless $access_key;
    
    ok $scheduler = VRPipe::Scheduler->create(type => 'ec2'), q[able to get the ec2 scheduler using get(type => 'ec2')];
    is $scheduler->type, 'ec2', 'the type really is ec2';
    
    $requirements = VRPipe::Requirements->create(memory => 100, time => 120);
    is $scheduler->determine_queue($requirements), 't1.micro', 'determine_queue() gave t1.micro instance for 100MB and 2mins';
    
    $requirements = VRPipe::Requirements->create(memory => 1800, time => 120);
    is $scheduler->determine_queue($requirements), 'm1.medium', 'determine_queue() gave m1.medium instance for 1800MB and 2mins';
    
    my $scheduler_cmd_line = $scheduler->submit_command(
        requirements => $requirements,
        stdo_file    => '/dev/null',
        stde_file    => '/dev/null',
        cmd          => 'the cmd to run'
    );
    my $expected = q{perl .+ "VRPipe::Scheduler->get\(type => q\[ec2\]\)->scheduler_instance->submit\(@ARGV\)" queue \S+ memory 1800 count 1 cmd 'the cmd to run'};
    like $scheduler_cmd_line, qr/$expected/, 'the expected scheduler cmd line could be constructed using submit_command()';
}

# sge
ok $scheduler = VRPipe::Scheduler->create(type => 'sge'), q[able to get the sge scheduler using get(type => 'sge')];
is $scheduler->type, 'sge', 'the type really is sge';
SKIP: {
    {
        no warnings "exec";
        my $qconf_out = `qconf -help`;
        skip "SGE 8 is not installed", 1 unless $qconf_out && $qconf_out =~ /^SGE 8/;
    }
    
    $requirements = VRPipe::Requirements->create(memory => 1800, time => 120);
    my $scheduler_cmd_line = $scheduler->submit_command(
        requirements => $requirements,
        stdo_file    => '/dev/null',
        stde_file    => '/dev/null',
        cmd          => 'the cmd to run',
        count        => 5
    );
    my $expected = q{qsub -N \S+ -o /dev/null -e /dev/null -m n -l (?:h_data|h_rss|h_vmem|mem_free|s_data|s_rss|s_vmem|virtual_free)=1800M (?:-l [hs]_rt=120 )?-t 1-5 -V -cwd -b yes the cmd to run};
    like $scheduler_cmd_line, qr/$expected/, 'the expected scheduler cmd line could be constructed using submit_command()';
}

done_testing;
exit;
