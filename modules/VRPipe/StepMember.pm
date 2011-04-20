use VRPipe::Base;

class VRPipe::StepMember extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'step' => (is => 'rw',
                   isa => Persistent,
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'pipeline' => (is => 'rw',
                       isa => Persistent,
                       coerce => 1,
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1);
    
    __PACKAGE__->make_persistent(belongs_to => [pipeline => 'VRPipe::Pipeline'],
                                 has_one => [step => 'VRPipe::Step']);
}

1;