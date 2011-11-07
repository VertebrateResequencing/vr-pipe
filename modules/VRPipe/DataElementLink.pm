use VRPipe::Base;

class VRPipe::DataElementLink extends VRPipe::Persistent {
    
    has 'parent' => (is => 'rw',
                     isa => Persistent,
                     coerce => 1,
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1,
                     belongs_to => 'VRPipe::DataElement');
    
    has 'child' => (is => 'rw',
                    isa => Persistent,
                    coerce => 1,
                    traits => ['VRPipe::Persistent::Attributes'],
                    is_key => 1,
                    belongs_to => 'VRPipe::DataElement');
    
    __PACKAGE__->make_persistent();
}

1;