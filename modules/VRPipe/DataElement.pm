use VRPipe::Base;

class VRPipe::DataElement extends VRPipe::Persistent {
    has 'datasource' => (is => 'rw',
                         isa => Persistent,
                         coerce => 1,
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_key => 1,
                         belongs_to => 'VRPipe::DataSource');
    
    has 'result' => (is => 'rw',
                     isa => 'HashRef',
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'withdrawn' => (is => 'rw',
                        isa => 'Bool',
                        traits => ['VRPipe::Persistent::Attributes'],
                        default => 0);
    
    __PACKAGE__->make_persistent(); # has_many => [element_states => 'VRPipe::DataElementState'] doesn't work because of ordering issues?
    
    method element_states {
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('DataElementState')->search({ dataelement => $self });
        my @states;
        while (my $state = $rs->next) {
            push(@states, $state);
        }
        return @states;
    }
    
    method start_from_scratch (VRPipe::PipelineSetup $setup, ArrayRef[PositiveInt] $step_numbers?) {
        VRPipe::DataElementState->get(pipelinesetup => $setup, dataelement => $self)->start_from_scratch($step_numbers ? $step_numbers : ());
    }
}

1;