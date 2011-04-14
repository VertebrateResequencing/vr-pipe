=head1 NAME

VRPipe::Base::SpawnProcesses - a role for sub-process management

=head1 SYNOPSIS

use VRPipe::Base;

class VRPipe::NeedsToSpawnAProcess with VRPipe::Base::SpawnProcesses {
    use MooseX::Workers::Job;

    method _build_processes {
        return MooseX::Workers::Job->new(name => $name,
                                         command => sub { #... });
    }
}


package main;

use VRPipe::NeedsToSpawnAProcess;

my $spawner = VRPipe::NeedsToSpawnAProcess->new();
$spawner->run;


=head1 DESCRIPTION

This is MooseX::Workers with a VRPipe-specific interface that allows you to
easily run one or more child processes.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use VRPipe::Base;

role VRPipe::Base::SpawnProcesses with MooseX::Workers {
    requires '_build_processes';
    
    has processes => (
        is      => 'ro',
        isa     => ArrayRefOfMWJobs,
        required => 1,
        builder => '_build_processes',
        coerce => 1
    );
    
    method run {
        foreach my $process (@{$self->processes}) {
            $self->run_command($process);
        }
        POE::Kernel->run();
    }
    
    method worker_manager_start { warn 'started worker manager' }
    method worker_manager_stop  { warn 'stopped worker manager'; }
    method max_workers_reached { warn 'maximum worker count reached' }
    method worker_error { warn "worker error: ", join ' ', @_;  }
    method worker_done (MooseX::Workers::Job $job) { printf ("%s(%s,%s) finished\n", $job->name, $job->ID, $job->PID); }
    method worker_started (MooseX::Workers::Job $job) { printf ("%s(%s,%s) started\n", $job->name, $job->ID, $job->PID); }
    
    method sig_TERM { warn 'Handled TERM' }
    method sig_child { warn "sig child: ", join ' ', @_; }
}

1;