use VRPipe::Base;

class VRPipe::StepStats extends VRPipe::Persistent {
    has 'step' => (is => 'rw',
                   isa => Persistent,
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1,
                   belongs_to => 'VRPipe::Step');
    
    has 'pipelinesetup' => (is => 'rw',
                            isa => Persistent,
                            coerce => 1,
                            traits => ['VRPipe::Persistent::Attributes'],
                            is_key => 1,
                            belongs_to => 'VRPipe::PipelineSetup');
    
    has 'submission' => (is => 'rw',
                         isa => Persistent,
                         coerce => 1,
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_key => 1,
                         belongs_to => 'VRPipe::Submission');
    
    has 'memory' => (is => 'rw',
                     isa => IntSQL[6],
                     traits => ['VRPipe::Persistent::Attributes']);
    
    has 'time' => (is => 'rw',
                   isa => IntSQL[6],
                   traits => ['VRPipe::Persistent::Attributes']);
    
    __PACKAGE__->make_persistent();
}

1;