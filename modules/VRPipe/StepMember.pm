use VRPipe::Base;

class VRPipe::StepMember extends VRPipe::Persistent {
    has 'step' => (is => 'rw',
                   isa => Persistent,
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1,
                   belongs_to => 'VRPipe::Step');
    
    has 'pipeline' => (is => 'rw',
                       isa => Persistent,
                       coerce => 1,
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       belongs_to => 'VRPipe::Pipeline');
    
    has 'step_number' => (is => 'rw',
                          isa => IntSQL[4],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1);
    
    __PACKAGE__->make_persistent();
    
    around step (Defined :$previous_step_outputs?, VRPipe::StepState :$step_state?) {
        my $step = $self->$orig();
        $step->step_state($step_state) if $step_state;
        $step->previous_step_outputs($previous_step_outputs) if $previous_step_outputs;
        return $step;
    }
}

1;