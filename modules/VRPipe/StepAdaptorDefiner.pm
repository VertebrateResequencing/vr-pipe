use VRPipe::Base;

class VRPipe::StepAdaptorDefiner {
    has 'from_step' => (is => 'ro',
                        isa => 'Int');
    
    has 'to_step' => (is => 'ro',
                      isa => PositiveInt,
                      coerce => 1);
    
    has 'from_key' => (is => 'ro',
                       isa => 'Str',
                       builder => '_from_key_builder',
                       lazy => 1);
    
    has 'to_key' => (is => 'ro',
                     isa => 'Str');
    
    method _from_key_builder {
        if ($self->from_step == 0) {
            return 'data_element';
        }
        else {
            $self->throw("from_key must be supplied if from_step is not 0");
        }
    }
    
    method define (Persistent|VRPipe::Pipeline $pipeline) {
        my $sa = VRPipe::StepAdaptor->get(pipeline => $pipeline, to_step => $self->to_step);
        my $adaptor_hash = $sa->adaptor_hash;
        $adaptor_hash->{$self->to_key}->{$self->from_key} = $self->from_step;
        $sa->adaptor_hash($adaptor_hash);
        $sa->update;
        return $sa;
    }
}

1;