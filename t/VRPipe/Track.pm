use VRPipe::Base;

class t::VRPipe::Track extends VRPipe::Persistent {
    has 'cd' => (is => 'rw',
                 isa => Persistent,
                 traits => ['VRPipe::Persistent::Attributes'],
                 belongs_to => 't::VRPipe::CD');
    
    has 'title' => (is => 'rw',
                    isa => Varchar[64],
                    traits => ['VRPipe::Persistent::Attributes']);
    
    __PACKAGE__->make_persistent();
}

1;
