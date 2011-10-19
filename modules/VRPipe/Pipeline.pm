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
    
    __PACKAGE__->make_persistent(has_many => [steps => 'VRPipe::StepMember']); # *** do not call steps...
    
    # for some bizarre reason, calling a method called steps (regardless of if
    # the has_many option is given to make_persistent) causes a memory leak,
    # so we define a method called step_members that does pretty much exactly
    # what steps method did
    sub step_members {
        my $self = shift;
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('StepMember')->search({ pipeline => $self->id });
        
        my @stepms;
        while (my $stepm = $rs->next) {
            push(@stepms, $stepm);
        }
        
        return @stepms;
    }
}

1;