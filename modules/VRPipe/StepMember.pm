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
                   is_key => 1,
                   belongs_to => 'VRPipe::Step');
    
    has 'pipeline' => (is => 'rw',
                       isa => Persistent,
                       coerce => 1,
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       belongs_to => 'VRPipe::Pipeline');
    
    __PACKAGE__->make_persistent();
    
    around step (Defined :$previous_step_outputs?, VRPipe::StepState :$step_state?) {
        my $step = $self->$orig();
        $step->step_state($step_state) if $step_state;
        $step->previous_step_outputs($previous_step_outputs) if $previous_step_outputs;
        return $step;
    }
}

1;