use VRPipe::Base;

class VRPipe::StepBehaviourDefiner {
    use Data::Compare;
    
    has 'after_step' => (is => 'ro',
                         isa => PositiveInt);
    
    has 'behaviour' => (is => 'ro',
                        isa => 'Str');
    
    has 'act_on_steps' => (is => 'ro',
                           isa => 'ArrayRef');
    
    has 'regulated_by' => (is => 'ro',
                           isa => 'Str');
    
    has 'default_regulation' => (is => 'ro',
                                 isa => 'Bool');
    
    method define (Persistent|VRPipe::Pipeline $pipeline) {
        my $sb = VRPipe::StepBehaviour->get(pipeline => $pipeline, after_step => $self->after_step, behaviour => $self->behaviour);
        my $array = $sb->behaviour_array;
        
        my $behaviour = $self->behaviour;
        my @steps = sort { $a <=> $b } @{$self->act_on_steps};
        my $already_have = 0;
        foreach my $existing (@$array) {
            my ($this_b, @these_steps) = @$existing;
            next unless $this_b eq $behaviour;
            
            if (Compare(\@steps, \@these_steps)) {
                $already_have = 1;
                last;
            }
        }
        
        unless ($already_have) {
            push(@$array, [$behaviour, @steps]);
        }
        
        $sb->behaviour_array($array);
        $sb->regulated_by($self->regulated_by);
        $sb->default_regulation($self->default_regulation);
        $sb->update;
        return $sb;
    }
}

1;