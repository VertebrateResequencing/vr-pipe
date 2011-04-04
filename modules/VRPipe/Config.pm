=pod

my $config = VRPipe::Config->new( config_file => '/path/to/config/myapp.ini' );
$config->write_file( ... );

=cut

use VRPipe::Base;

class VRPipe::Config {
    use VRPipe::Base::Configuration;
    
    my $question_number = 0;
    
    has production_dbtype => (
        is      => 'rw',
        default => 'mysql',
        question => 'What DRM should be used for production?',
        valid => [qw(mysql postgres)],
        question_number => ++$question_number
    );
    
    has production_dbname => (
        is      => 'rw',
        default => 'db_name_default_from_Config',
        question => 'What is the name of your production database?',
        question_number => ++$question_number
    );
    
    has production_username => (
        is      => 'rw',
        default => 'db_user_default_from_Config',
        documentation => 'What username is used to connect to your production database?',
        env     => 'VRTRACK_RW_USER',
        question_number => ++$question_number
    );
    
    has testing_dbtype => (
        is      => 'rw',
        default => 'sqlite',
        question => 'What DRM should be used for testing?',
        valid => [qw(sqlite mysql postgres)],
        question_number => ++$question_number
    );
    
    has testing_dbname => (
        is      => 'rw',
        default => ':memory:',
        question => 'What is the name of your testing database?',
        question_number => ++$question_number
    );
    
    has testing_username => (
        is      => 'rw',
        default => '',
        documentation => 'What username is used to connect to your testing database?',
        question_number => ++$question_number
    );
}

1;