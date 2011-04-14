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
        return MooseX::Workers::Job->new(
            name    => 'heartbeat',
            command => sub { while (1) { sleep 2; print "heartbeat\n"; } }
        );
    }
}

package main;

use VRPipe::NeedsToSpawnAProcess;

my $spawner = VRPipe::NeedsToSpawnAProcess->new();
$spawner->run;

=head1 DESCRIPTION

This role sets up an initial heart-beat process that will run prior to the
normal processes and gets terminated when all others are done.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use VRPipe::Base;

role VRPipe::Base::LivingProcesses with VRPipe::Base::SpawnProcesses {
    requires '_build_heartbeat_sub';
    
    has heartbeat => (
        is      => 'ro',
        isa     => 'MooseX::Workers::Job',
        builder => '_build_heartbeat_sub',
        lazy    => 1,
        required => 1
    );
    
    has heartbeat_interval => (
        is      => 'ro',
        isa     => 'Int',
        default => 60
    );
    
    before run {
        $self->run_command($self->heartbeat);
    }
    
    after sig_child {
        # kill the heartbeat
        if ($self->num_workers == 1) {
            #$self->kill_worker(??); no idea what the args to kill_worker are supposed to be!
            $self->Engine->get_worker(1)->kill();
        }
    }
}

1;