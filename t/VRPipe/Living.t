#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use EV;
use AnyEvent;

BEGIN {
    use Test::Most tests => 22;
    use VRPipeTest;
    use TestPipelines;
}

# Test core-functionality common to all Living objects
{
    # we use Runner as the example, which is exactly like the Living base
    # except that it has some extra methods
    my $l = VRPipe::Runner->create(cmd => 'a');
    $l->start_beating;
    my ($heartbeat_interval, $survival_time, $sleep);
    
    # basic heartbeat tests
    my @watchers;
    $watchers[0] = EV::timer 0, 0, sub {
        is $l->alive, 1, 'created a Living object and it is alive';
        my $l_count = VRPipe::Runner->search({});
        is $l_count, 1, 'we have an entry in the db';
        ok my $heartbeat_time = $l->heartbeat->epoch, 'it had a heartbeat';
        ok $heartbeat_interval = $l->heartbeat_interval, 'it had a heartbeat_interval';
        ok $survival_time      = $l->survival_time,      'it had a survival_time';
        cmp_ok $survival_time, '>', $heartbeat_interval, 'survival time is greater than the heartbeat_interval';
        
        $sleep = $heartbeat_interval + 2;
        $watchers[1] = EV::timer $sleep, 0, sub {
            my $new_heartbeat_time = $l->heartbeat->epoch;
            cmp_ok $new_heartbeat_time, '>', $heartbeat_time, 'the heart is beating';
            $l->stop_beating;
            
            $watchers[2] = EV::timer $sleep, 0, sub {
                cmp_ok $l->heartbeat->epoch, '==', $new_heartbeat_time, 'after calling stop_beating, our heart stopped beating';
                is $l->alive, 1, 'but right now we are still considered alive';
                
                my $time_until_dead = $survival_time - $sleep + 2;
                $watchers[3] = EV::timer $time_until_dead, 0, sub {
                    is $l->alive, 0, 'after waiting longer than survival_time without a heartbeat, we died';
                    $l_count = VRPipe::Runner->search({});
                    is $l_count, 0, 'our entry was removed from the db';
                    EV::unloop;
                };
            };
        };
    };
    EV::run;
    @watchers = ();
    
    # test that we can see if a Living created in another process is still
    # alive, and what happens to it if we suspend that process
    my $tmp_dir           = $l->tempdir;
    my $child_stderr_file = file($tmp_dir, 'child_stderr');
    my $child_pid         = fork();
    if (!defined $child_pid) {
        die "attempt to fork failed: $!";
    }
    elsif ($child_pid) {
        # parent
        $watchers[0] = EV::timer $sleep, 0, sub {
            my $l_count = VRPipe::Runner->search({ cmd => 'b' });
            is $l_count, 1, 'child process created a Living object and parent was able to find it';
            my $test_l = VRPipe::Runner->get(cmd => 'b');
            is $test_l->alive, 1, 'it is alive';
            
            # suspend the child for longer than survival_time
            kill('SIGSTOP', $child_pid);
            my $time_until_dead = $survival_time + 2;
            $watchers[1] = EV::timer $time_until_dead, 0, sub {
                is $test_l->alive, 0, 'it died after a while of being suspended';
                
                # now resume the child and see what happens...
                kill('SIGCONT', $child_pid);
                EV::unloop;
            };
        };
        EV::run;
    }
    elsif ($child_pid == 0) {
        # child
        open(STDERR, '>', $child_stderr_file);
        my $child_l = VRPipe::Runner->create(cmd => 'b');
        $child_l->start_beating;
        
        my $timeout = $sleep + $survival_time + 4;
        $watchers[2] = EV::timer $timeout, 0, sub {
            warn "child failed to detect being murdered; we timed out instead\n";
            EV::unloop;
        };
        
        EV::run;
        
        exit(0);
    }
    waitpid($child_pid, 0);
    @watchers = ();
    my $child_stderr = file($child_stderr_file)->slurp;
    like $child_stderr, qr/We were murdered by another process/, 'child detects parent murdered it and ended';
    unlike $child_stderr, qr/we timed out instead/, 'child did not exit due to the timeout';
    undef $l;
}

# FarmServer-specific tests
{
    # make some pipelinesetups we can claim
    my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_stats', 'bam_files');
    my $ds = VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data datasource.bam_fofn))->absolute);
    my %common_setup_args = (datasource => $ds, output_root => $output_dir, pipeline => $pipeline, options => {});
    VRPipe::PipelineSetup->create(name => 'one',   %common_setup_args);
    VRPipe::PipelineSetup->create(name => 'two',   %common_setup_args);
    VRPipe::PipelineSetup->create(name => 'three', %common_setup_args, desired_farm => 'bar');
    VRPipe::PipelineSetup->create(name => 'four',  %common_setup_args, desired_farm => 'foo');
    VRPipe::PipelineSetup->create(name => 'five',  %common_setup_args, desired_farm => 'foo');
    
    is_deeply [map { $_->controlling_farm } VRPipe::PipelineSetup->search({})], [undef, undef, undef, undef, undef], 'to start with, we have 5 uncontrolled setups';
    
    ok my $fs = VRPipe::FarmServer->create(farm => 'foo', only_ours => 1), 'able to create a FarmServer';
    my $num_setups = $fs->claim_setups;
    is $num_setups, 2, 'claim_setups returned the correct number of setups when only_ours was used';
    $fs->commit_suicide(no_die => 1);
    undef $fs;
    my $fs_count = VRPipe::FarmServer->search({});
    is $fs_count, 0, 'commit_suicide worked, removing the farmserver from the db';
    $fs = VRPipe::FarmServer->create(farm => 'foo');
    $num_setups = $fs->claim_setups;
    is $num_setups, 4, 'claim_setups returned the correct number of setups when only_ours was not used';
    undef $fs;
}

# Runner-specific tests
{
    ok my $r = VRPipe::Runner->create(cmd => 'ls'), 'able to create a Runner';
}

exit;
