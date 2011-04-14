use VRPipe::Base;

class t::VRPipe::SimpleSpawner with VRPipe::Base::SpawnProcesses {
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
        return MooseX::Workers::Job->new(name => 'simple',
                                         command => sub { print POE::Filter::Reference->new->put([ {msg => "simple"} ]);
                                                          print "foo\n"; });
    }
    
    method stdout_filter { new POE::Filter::Reference }
    #
    #method worker_stdout (HashRef $output) {
    #    warn "calling add_output($output)\n";
    #    $self->add_output($output->{msg});
    #    warn "output is now ", join(",", @{$self->output}), "\n";
    #}
    
    method worker_stdout (Str $output, MooseX::Workers::Job $job) {
        print $job->name, "(", $job->ID, ",", $job->PID, ") said ", $output, "\n";
    }
}

1;