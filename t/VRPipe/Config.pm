use VRPipe::Base;

class t::VRPipe::Config {
    use VRPipe::Base::Configuration;
    
    my $question_number = 0;
    
    has first_key => (
        is      => 'rw',
        default => 'foo',
        question => 'First question?',
        valid => [qw(foo bar)],
        question_number => ++$question_number
    );
    
    has second_key => (
        is      => 'rw',
        question => 'Second question?',
        question_number => ++$question_number
    );
    
    has third_key => (
        is      => 'rw',
        env     => 'TVRPIPE_THIRDKEY',
        question_number => ++$question_number
    );
}

1;