use VRPipe::Base;

class VRPipe::DataElement extends VRPipe::Persistent {
    has 'datasource' => (is => 'rw',
                         isa => Persistent,
                         coerce => 1,
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_key => 1,
                         belongs_to => 'VRPipe::DataSource');
    
    has 'result' => (is => 'rw',
                     isa => Varchar[256],
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    __PACKAGE__->make_persistent();
}

1;