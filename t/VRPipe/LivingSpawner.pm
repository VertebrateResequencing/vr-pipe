use VRPipe::Base;

class t::VRPipe::LivingSpawner with VRPipe::Base::LivingProcesses {
    use MooseX::Workers::Job;
    
    method _build_processes {
        return MooseX::Workers::Job->new(name => 'living',
                                         command => sub { sleep 1; print "living\n"; });
    }
    
    method _build_heartbeat_sub {
        return sub { print "heartbeat\n" }
    }
    
    method worker_stdout (Str $output, MooseX::Workers::Job $job) {
        print $job->name, " said ", $output, "\n";
    }
}

1;