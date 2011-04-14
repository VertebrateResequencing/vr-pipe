use VRPipe::Base;

class t::VRPipe::LivingSpawner with VRPipe::Base::LivingProcesses {
    use MooseX::Workers::Job;
    use POE qw(Filter::Reference);
    
    has output => (
        traits  => ['Array'],
        is      => 'ro',
        isa     => 'ArrayRef[Str]',
        default => sub { [] },
        handles => {
            all_outputs    => 'elements',
            add_output     => 'push',
            map_outputs    => 'map',
            output_count   => 'count',
            sorted_outputs => 'sort',
        },
    );
    
    method _build_processes {
        return MooseX::Workers::Job->new(name => 'living',
                                         command => sub { print POE::Filter::Reference->new->put([ {msg => "living"} ]); });
    }
    
    method _build_heartbeat_sub {
        return MooseX::Workers::Job->new(
            name    => 'heartbeat',
            command => sub { while (1) { sleep 2; print POE::Filter::Reference->new->put([ {msg => "heartbeat"} ]); } }
        );
    }
    
    sub stdout_filter { new POE::Filter::Reference }
    
    method worker_stdout (HashRef $output) {
        $self->add_output($output->{msg});
    }
}

1;