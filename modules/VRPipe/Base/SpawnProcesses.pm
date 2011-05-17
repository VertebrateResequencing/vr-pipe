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
    use Sys::CPU;
    my $cpu_count = Sys::CPU::cpu_count();
    if ($cpu_count < 2) {
        $cpu_count = 2;
    }
    
    requires '_build_processes';
    
    has processes => (
        is      => 'ro',
        isa     => ArrayRefOfMWJobs,
        required => 1,
        builder => '_build_processes',
        coerce => 1
    );
    
    has max_processes => (
        is      => 'rw',
        isa     => 'Int',
        default => $cpu_count,
    );
    
    sub catch_run_interupts {
        my $signame = shift;
        die "Somebody sent me a SIG$signame";
    }
    
    method run {
        local $SIG{INT} = \&catch_run_interupts;
        local $SIG{TERM} = \&catch_run_interupts;
        
        $self->max_workers($self->max_processes);
        foreach my $process (@{$self->processes}) {
            $self->enqueue($process);
        }
        POE::Kernel->run();
    }
    
    # so that we can use method modifiers in "inheriting" roles and classes
    method worker_manager_start {  }
    method worker_manager_stop  {  }
    method max_workers_reached {  }
    method worker_error {  }
    method worker_started (MooseX::Workers::Job $job) {
        #warn sprintf("%s(%s,%s) started\n", $job->name, $job->ID, $job->PID);
    }
    method worker_done (MooseX::Workers::Job $job) {
        #warn sprintf("%s(%s,%s) finished\n", $job->name, $job->ID, $job->PID);
    }
    method sig_TERM {  }
    method sig_child {  }
}

1;