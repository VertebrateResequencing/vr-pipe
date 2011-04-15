=head1 NAME

VRPipe::Base::LivingProcesses - a role for sub-processes with heartbeats

=head1 SYNOPSIS

use VRPipe::Base;

class VRPipe::NeedsToSpawnAProcess with VRPipe::Base::LivingProcesses {
    use MooseX::Workers::Job;

    method _build_processes {
        return MooseX::Workers::Job->new(name => $name,
                                         command => sub { #... });
    }
    
    method _build_heartbeat_sub {
        return sub { print "heartbeat\n"; }
    }
}

package main;

use VRPipe::NeedsToSpawnAProcess;

my $spawner = VRPipe::NeedsToSpawnAProcess->new();
$spawner->run;

=head1 DESCRIPTION

This role sets up an initial heart-beat process that will run prior to the
normal processes and that gets terminated when all others are done.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use VRPipe::Base;

role VRPipe::Base::LivingProcesses with VRPipe::Base::SpawnProcesses {
    requires '_build_heartbeat_sub';
    has heartbeat => (
        is      => 'ro',
        isa     => 'CodeRef',
        builder => '_build_heartbeat_sub',
        lazy    => 1,
        required => 1
    );
    
    has heartbeat_interval => (
        is      => 'ro',
        isa     => 'Int',
        default => 60
    );
    
    has _heartbeat_mwj => (
        is      => 'ro',
        isa     => 'MooseX::Workers::Job',
        builder => '_build_heartbeat_mwj',
        lazy    => 1
    );
    
    method _build_heartbeat_mwj {
        my $heartbeat_sub = $self->heartbeat;
        my $interval = $self->heartbeat_interval;
        return MooseX::Workers::Job->new(
            name    => 'heartbeat',
            command => sub { while (1) { sleep $interval; &$heartbeat_sub; } }
        );
    }
    
    before run {
        if ($self->max_processes < 2) {
            $self->max_processes(2);
            $self->max_workers(2);
        }
        $self->enqueue($self->_heartbeat_mwj);
    }
    
    after sig_child {
        # kill the heartbeat
        if ($self->num_workers == 1) {
            my $engine = $self->Engine;
            my @ids = sort { $a <=> $b } $engine->get_worker_ids;
            $engine->get_worker($ids[0])->kill();
        }
    }
}

1;