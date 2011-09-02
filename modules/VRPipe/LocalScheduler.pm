use VRPipe::Base;

class VRPipe::LocalScheduler {
    with 'MooseX::Daemonize';
    use VRPipe::LocalSchedulerJob;
    
    after start {
        return unless $self->is_daemon;
        # your daemon code here ...
    }
    
    method submit {
        warn "submit received\n";
    }
}