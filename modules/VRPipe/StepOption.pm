use VRPipe::Base;

class VRPipe::StepOption extends VRPipe::Persistent {
    has 'description' => (is => 'rw',
                          isa => Varchar[256],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1);
    
    has 'optional' => (is => 'rw',
                       isa => 'Bool',
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       default => 0,
                       allow_key_to_default => 1);
    
    has 'default_value' => (is => 'rw',
                            isa => Varchar[256],
                            traits => ['VRPipe::Persistent::Attributes'],
                            is_key => 1,
                            default => '',
                            allow_key_to_default => 1);
    
    has 'allowed_values' => (is => 'rw',
                            isa => 'ArrayRef',
                            traits => ['VRPipe::Persistent::Attributes'],
                            is_key => 1,
                            default => sub { [] },
                            allow_key_to_default => 1);
    
    __PACKAGE__->make_persistent();
}

1;