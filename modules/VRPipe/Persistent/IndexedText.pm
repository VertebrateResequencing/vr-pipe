use VRPipe::Base;

class VRPipe::Persistent::IndexedText extends VRPipe::Persistent {
    has 'sha' => (is => 'rw',
                  isa => Char[43],
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1);
    
    has 'text' => (is => 'rw',
                   isa => Text,
                   traits => ['VRPipe::Persistent::Attributes']);
    
    __PACKAGE__->make_persistent();
}

1;