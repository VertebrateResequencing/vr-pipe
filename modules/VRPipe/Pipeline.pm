use VRPipe::Base;

class VRPipe::Pipeline extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'description' => (is => 'rw',
                         isa => Varchar[256],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    __PACKAGE__->make_persistent(has_many => [steps => 'VRPipe::StepMember']);
}

1;