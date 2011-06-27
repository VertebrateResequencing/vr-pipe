use VRPipe::Base;

class t::VRPipe::Artist extends VRPipe::Persistent {
    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'age' => (is => 'rw',
                   isa => IntSQL[3],
                   default => 99,
                   traits => ['VRPipe::Persistent::Attributes']);
    
    __PACKAGE__->make_persistent(has_many => [cds => 't::VRPipe::CD']);
}

1;
