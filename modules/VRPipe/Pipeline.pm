use VRPipe::Base;

class VRPipe::Pipeline extends VRPipe::Persistent with VRPipe::PipelineRole {
    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has '_num_steps' => (is => 'rw',
                         isa => IntSQL[4],
                         traits => ['VRPipe::Persistent::Attributes'],
                         default => 0);
    
    has 'description' => (is => 'rw',
                          isa => Text,
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1);
    
    __PACKAGE__->make_persistent(has_many => [steps => 'VRPipe::StepMember']);
}

1;