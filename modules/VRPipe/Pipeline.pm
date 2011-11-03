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
    
    # steps must be called to initially create stepmembers for a new pipeline,
    # but currently causes a memory leak, so we want to call it only once.
    # Everywhere else we call step_members to just get back the existing
    # step members
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