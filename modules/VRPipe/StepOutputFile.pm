use VRPipe::Base;

class VRPipe::StepOutputFile extends VRPipe::Persistent {
    has 'stepstate' => (is => 'rw',
                        isa => Persistent,
                        coerce => 1,
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_key => 1,
                        belongs_to => 'VRPipe::StepState');
    
    has 'file' => (is => 'rw',
                   isa => Persistent,
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1,
                   belongs_to => 'VRPipe::File');
    
    has 'output_key' => (is => 'rw',
                         isa => Varchar[64],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_key => 1);
    
    __PACKAGE__->make_persistent();
}

1;