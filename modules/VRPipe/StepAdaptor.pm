use VRPipe::Base;

class VRPipe::StepAdaptor extends VRPipe::Persistent {
    has 'pipeline' => (is => 'rw',
                       isa => Persistent,
                       coerce => 1,
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       belongs_to => 'VRPipe::Pipeline');
    
    has 'to_step' => (is => 'rw',
                      isa => IntSQL[4],
                      traits => ['VRPipe::Persistent::Attributes'],
                      is_key => 1);
    
    has 'adaptor_hash' => (is => 'rw',
                           isa => 'HashRef',
                           traits => ['VRPipe::Persistent::Attributes'],
                           default => sub { {} });
    
    __PACKAGE__->make_persistent();
    
    method adapt (Str :$input_key!, PreviousStepOutput :$pso?, VRPipe::DataElement :$data_element?) {
        if ((defined $pso ? 1 : 0) + (defined $data_element ? 1 : 0) != 1) {
            $self->throw("Exactly one of pso or data_element must be supplied");
        }
        
        my $hash = $self->adaptor_hash;
        
        if (defined $hash->{$input_key}) {
            my $ref = $hash->{$input_key};
            
            my @results;
            while (my ($key, $from_step) = each %$ref) {
                if ($data_element) {
                    if ($key eq 'data_element' && $from_step == 0) {
                        my $result = $data_element->result;
                        my $paths = $result->{paths} || $self->throw("data element ".$data_element->id." gave a result with no paths");
                        foreach my $path (@$paths) {
                            push(@results, VRPipe::File->get(path => file($path)->absolute));
                        }
                    }
                }
                else {
                    my $result = $pso->{$key}->{$from_step} || next;
                    push(@results, (ref($result) && ref($result) eq 'ARRAY') ? @$result : $result);
                }
            }
            
            return @results ? \@results : undef;
        }
        
        return;
    }
}

1;