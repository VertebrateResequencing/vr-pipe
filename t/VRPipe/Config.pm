use VRPipe::Base;

class t::VRPipe::Config {
    use VRPipe::Base::Configuration;
    
    my $question_number = 0;
    
    has first_key => (
        is              => 'rw',
        default         => 'foo',
        question        => 'First question?',
        valid           => [qw(foo bar make_second_skip)],
        question_number => ++$question_number
    );
    
    has second_key => (
        is              => 'rw',
        question        => 'Second question?',
        skip            => '_skip_second_key',
        question_number => ++$question_number
    );
    
    has third_key => (
        is              => 'rw',
        env             => 'TVRPIPE_THIRDKEY',
        builder         => '_build_third_key',
        question_number => ++$question_number
    );
    
    has secret_key => (
        is              => 'rw',
        secure          => 1,
        question_number => ++$question_number
    );
    
    method _skip_second_key {
        my $first_key = $self->first_key;
        if ($first_key && $first_key eq 'make_second_skip') {
            return 1;
        }
        return 0;
    }
    
    method _build_third_key {
        if ($self->_skip_second_key) {
            return 'based on skipped second key';
        }
        return;
    }
}

1;
