use VRPipe::Base;

class VRPipe::Step extends VRPipe::Persistent with VRPipe::StepRole {
    use VRPipe::StepNonPersistentFactory;
    
    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'inputs_definition' => (is => 'rw',
                                isa => PersistentHashRef,
                                traits => ['VRPipe::Persistent::Attributes']);
    
    has 'body_sub' => (is => 'rw',
                       isa => 'CodeRef',
                       traits => ['VRPipe::Persistent::Attributes']);
    
    has 'post_process_sub' => (is => 'rw',
                               isa => 'CodeRef',
                               traits => ['VRPipe::Persistent::Attributes']);
    
    has 'outputs_definition' => (is => 'rw',
                                 isa => PersistentHashRef,
                                 traits => ['VRPipe::Persistent::Attributes']);
    
    has 'description' => (is => 'rw',
                          isa => Varchar[64],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1);
    
    __PACKAGE__->make_persistent();
    
    around get (ClassName|Object $self: %args) {
        # first see if there is a corresponding VRPipe::Steps::$name class
        if (defined $args{name}) {
            try {
                eval "require VRPipe::Steps::$args{name};"; # eval within a try because can't get require to work with a variable name otherwise?!
                die "$@\n" if $@;
                my $obj = VRPipe::StepNonPersistentFactory->create($args{name}, {});
                return $obj;
            }
            catch ($err) {
                unless ($err =~ /^Can't locate/) {
                    $self->throw($err);
                }
            }
        }
        
        # otherwise, assume we're getting/creating a step from persistent db
        return $self->$orig(%args);
    }
}

1;