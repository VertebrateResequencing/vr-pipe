#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 17;
    #use Test::Output; complains about inability to redirect into tied STDOUT
    
    use_ok('t::VRPipe::SimpleSpawner');
    use_ok('t::VRPipe::LivingSpawner');
    use_ok('VRPipe::Cmd');
    use_ok('MooseX::Workers::Job');
}

# test the simple spawner, which uses the basic SpawnProcesses
ok my $ss = t::VRPipe::SimpleSpawner->new(), 'made a new SimpleSpawner using default processe';
is run($ss), "simple said simple\n", 'output from a single process was ok';

my @jobs;
foreach my $name ('foo', 'bar', 'henry') {
    my $job = MooseX::Workers::Job->new(
        name    => $name,
        command => sub { for (1..3) { sleep 1; print "$name\n"; } }
    );
    push(@jobs, $job);
}

ok $ss = t::VRPipe::SimpleSpawner->new(processes => \@jobs, max_processes => 4), 'made a new SimpleSpawner with custom processes';
is run($ss), join("\n", map { "$_ said $_" } qw(foo bar henry foo bar henry foo bar henry))."\n", 'stdout from 3 processes interleaved as expected';
is_deeply $ss->output, [qw(foo bar henry foo bar henry foo bar henry)], 'output from 3 processes could also be collated in the parent array ref';

# test the living spawner, which uses LivingProcesses
ok my $ls = t::VRPipe::LivingSpawner->new(heartbeat_interval => 1, max_processes => 4), 'made a new LivingSpawner using default processes';
is run($ls), "heartbeat said heartbeat\nliving said living\nheartbeat said heartbeat\n", 'output from a single process was ok';
ok $ls = t::VRPipe::LivingSpawner->new(processes => \@jobs, heartbeat_interval => 1, max_processes => 4), 'made a new LivingSpawner with custom processes';
is run($ls), join("\n", map { "$_ said $_" } qw(heartbeat foo bar henry heartbeat foo bar henry heartbeat foo bar henry heartbeat))."\n", 'stdout from 3 processes and the heartbeat interleaved as expected';

# test behaviour when max_workers is set low
$ss = t::VRPipe::SimpleSpawner->new(processes => [$jobs[0], $jobs[1]]);
cmp_ok $ss->max_processes, '>=', 2, 'default max_processes is 2 or more';
$ss = t::VRPipe::SimpleSpawner->new(processes => [$jobs[0], $jobs[1]], max_processes => 1);
is $ss->max_processes, 1, 'max_processes could be set to 1';
is run($ss), join("\n", map { "$_ said $_" } qw(foo foo foo bar bar bar))."\n", 'stdout from 2 processes are sequential as expected when only 1 process runs at a time';
$ls = t::VRPipe::LivingSpawner->new(processes => [$jobs[0], $jobs[1]], heartbeat_interval => 1, max_processes => 1);
is run($ls), join("\n", map { "$_ said $_" } qw(heartbeat foo heartbeat foo heartbeat foo
                                                heartbeat bar heartbeat bar heartbeat bar heartbeat))."\n", 'stdout from 2 processes and a heartbeat interleaved as expected when only 1 process runs at a time';


exit;

sub run {
    my $obj = shift;
    my $stdout = '';
    {
        local *STDOUT;
        open STDOUT, '>', \$stdout or die "Cannot open STDOUT to a scalar: $!";
        $obj->run;
        close STDOUT or die "Cannot close redirected STDOUT: $!";
    }
    return $stdout;
}