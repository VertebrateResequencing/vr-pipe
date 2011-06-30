use VRPipe::Base;

class VRPipe::StepAdaptor extends VRPipe::Persistent {
    has 'pipeline' => (is => 'rw',
                       isa => Persistent,
                       coerce => 1,
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       belongs_to => 'VRPipe::Pipeline');
    
    has 'before_step_number' => (is => 'rw',
                                 isa => IntSQL[4],
                                 traits => ['VRPipe::Persistent::Attributes'],
                                 is_key => 1);
    
    has 'adaptor_hash' => (is => 'rw',
                           isa => 'HashRef',
                           traits => ['VRPipe::Persistent::Attributes'],
                           default => sub { {} });
    
    __PACKAGE__->make_persistent();
    
    method adapt (Str :$input_key, PersistentFileHashRef :$pso?, VRPipe::DataElement :$data_element?) {
        if ((defined $pso ? 1 : 0) + (defined $data_element ? 1 : 0) != 1) {
            $self->throw("Exactly one of pso or data_element must be supplied");
        }
        
        my $hash = $self->adaptor_hash;
        if (defined $hash->{$input_key}) {
            my $allowed = $hash->{$input_key};
            my @allowed = ref($allowed) ? @$allowed : ($allowed);
            
            my @results;
            foreach my $key (@allowed) {
                if ($data_element) {
                    if ($key eq 'data_element') {
                        push(@results, VRPipe::File->get(path => file($data_element->result)->absolute));
                    }
                }
                else {
                    my $result = $pso->{$key} || next;
                    push(@results, (ref($result) && ref($result) eq 'ARRAY') ? @$result : $result);
                }
            }
            
            return @results ? \@results : undef;
        }
        
        return;
    }
}

1;