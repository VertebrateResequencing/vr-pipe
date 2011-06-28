use VRPipe::Base;

class VRPipe::Pipeline extends VRPipe::Persistent {
    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has '_num_steps' => (is => 'rw',
                        isa => IntSQL[4],
                        traits => ['VRPipe::Persistent::Attributes'],
                        default => 0);
    
    has 'description' => (is => 'rw',
                         isa => Varchar[256],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    method num_steps {
        return $self->_num_steps;
    }
    method _increment_steps (PositiveInt $step_num) {
        return $self->_num_steps($step_num);
    }
    
    __PACKAGE__->make_persistent(has_many => [steps => 'VRPipe::StepMember']);
    
    method add_step (VRPipe::Step $step) {
        my $step_num = $self->num_steps() + 1;
        VRPipe::StepMember->get(step => $step, pipeline => $self, step_number => $step_num);
        $self->_increment_steps($step_num);
        $self->update;
    }
}

1;