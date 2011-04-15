use VRPipe::Base;

class t::VRPipe::SimpleSpawner with VRPipe::Base::SpawnProcesses {
    use MooseX::Workers::Job;
    
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
                                         command => sub { print "simple\n"; });
    }
    
    method worker_stdout (Str $output, MooseX::Workers::Job $job) {
        print $job->name, " said ", $output, "\n";
        chomp($output);
        $self->add_output($output);
    }
}

1;