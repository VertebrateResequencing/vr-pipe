use VRPipe::Base;

class VRPipe::PersistentArrayMember extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'persistentarray' => (is => 'rw',
                              isa => Persistent,
                              coerce => 1,
                              traits => ['VRPipe::Persistent::Attributes'],
                              is_key => 1,
                              belongs_to => 'VRPipe::PersistentArray');
    
    has 'class' => (is => 'rw',
                    isa => Varchar[64],
                    traits => ['VRPipe::Persistent::Attributes'],
                    is_key => 1);
    
    has 'class_id' => (is => 'rw',
                       isa => Persistent, # uncoerced, we want this to be a plain int always
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1);
    
    __PACKAGE__->make_persistent();
}

1;