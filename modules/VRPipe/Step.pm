use VRPipe::Base;

class VRPipe::Step extends VRPipe::Persistent with VRPipe::StepRole {
    use VRPipe::StepNonPersistentFactory;
    
    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'options_definition' => (is => 'rw',
                                 isa => PersistentHashRef,
                                 traits => ['VRPipe::Persistent::Attributes'],
                                 default => sub { {} });
    
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
                          isa => Text,
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1);
    
    __PACKAGE__->make_persistent();
}

1;