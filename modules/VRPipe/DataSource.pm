use VRPipe::Base;

class VRPipe::DataSource extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'module' => (is => 'rw',
                     isa => Varchar[64],
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'method' => (is => 'rw',
                     isa => Varchar[64],
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'source' => (is => 'rw',
                     isa => Varchar[64],
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'options' => (is => 'rw',
                      isa => Varchar[64],
                      traits => ['VRPipe::Persistent::Attributes']);
    
    has 'description' => (is => 'rw',
                         isa => Varchar[64],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    __PACKAGE__->make_persistent(#has_many => [elements => 'VRPipe::DataElement']
                                 );
}

1;