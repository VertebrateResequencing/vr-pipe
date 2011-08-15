use VRPipe::Base;

class VRPipe::StepBehaviour extends VRPipe::Persistent {
    has 'pipeline' => (is => 'rw',
                       isa => Persistent,
                       coerce => 1,
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       belongs_to => 'VRPipe::Pipeline');
    
    has 'after_step' => (is => 'rw',
                         isa => IntSQL[4],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_key => 1);
    
    has 'behaviour_array' => (is => 'rw',
                              isa => 'ArrayRef',
                              traits => ['VRPipe::Persistent::Attributes'],
                              default => sub { [] });
    
    has 'regulated_by' => (is => 'rw',
                           isa => Varchar[64],
                           traits => ['VRPipe::Persistent::Attributes'],
                           default => '');
    
    has 'default_regulation' => (is => 'rw',
                                 isa => 'Bool',
                                 traits => ['VRPipe::Persistent::Attributes'],
                                 default => 0);
    
    has 'dataelement' => (is => 'rw',
                          isa => 'VRPipe::DataElement');
    has 'pipelinesetup' => (is => 'rw',
                            isa => 'VRPipe::PipelineSetup');
    
    __PACKAGE__->make_persistent();
    
    method behave (VRPipe::DataElement :$data_element, VRPipe::PipelineSetup :$pipeline_setup) {
        my $go_ahead = 0;
        my $regulated_by = $self->regulated_by;
        if ($regulated_by) {
            my $options = $pipeline_setup->options;
            if (defined $options->{$regulated_by}) {
                $go_ahead = $options->{$regulated_by};
            }
            else {
                $go_ahead = $self->default_regulation;
            }
        }
        else {
            $go_ahead = $self->default_regulation;
        }
        
        return unless $go_ahead;
        
        $self->dataelement($data_element);
        $self->pipelinesetup($pipeline_setup);
        
        my $array = $self->behaviour_array;
        
        foreach my $behaviour (@$array) {
            my ($method, @steps) = @$behaviour;
            $self->throw("'$method' is not a valid behaviour") unless $self->can($method);
            
            return $self->$method($self->step_numbers_to_states(\@steps));
        }
    }
    
    method step_numbers_to_states (ArrayRef[PositiveInt] $steps) {
        my $pipeline = $self->pipeline;
        my $de = $self->dataelement;
        my $ps = $self->pipelinesetup;
        
        my %step_nums = map { $_ => 1 } @$steps;
        
        my @states;
        foreach my $stepm ($pipeline->steps) {
            if (exists $step_nums{$stepm->step_number}) {
                push(@states, VRPipe::StepState->get(stepmember => $stepm, pipelinesetup => $ps, dataelement => $de));
            }
        }
        
        return \@states;
    }
    
    method delete_outputs (ArrayRef[VRPipe::StepState] $states) {
        foreach my $state (@$states) {
            $state->unlink_output_files;
        }
    }
    
    method start_over (ArrayRef[VRPipe::StepState] $states) {
        foreach my $state (@$states) {
            $state->start_over;
        }
    }
}

1;