use VRPipe::Base;

class VRPipe::Requirements extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'memory' => (is => 'rw',
                     isa => IntSQL[5],
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'time' => (is => 'rw',
                   isa => IntSQL[3],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'cpus' => (is => 'rw',
                   isa => IntSQL[5],
                   traits => ['VRPipe::Persistent::Attributes'],
                   default => 1,
                   is_key => 1,
                   allow_key_to_default => 1);
    
    has 'tmp_space' => (is => 'rw',
                        isa => IntSQL[5],
                        traits => ['VRPipe::Persistent::Attributes'],
                        default => 0,
                        is_key => 1,
                        allow_key_to_default => 1);
    
    has 'local_space' => (is => 'rw',
                          isa => IntSQL[5],
                          traits => ['VRPipe::Persistent::Attributes'],
                          default => 0,
                          is_key => 1,
                          allow_key_to_default => 1);
    
    has 'custom' => (is => 'rw',
                     isa => 'HashRef',
                     traits => ['VRPipe::Persistent::Attributes'],
                     default => sub { {} },
                     is_key => 1,
                     allow_key_to_default => 1);
    
    __PACKAGE__->make_persistent();
}

1;