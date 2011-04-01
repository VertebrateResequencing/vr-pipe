=pod

my $config = VRPipe::Config->new( config_file => '/path/to/config/myapp.ini' );
$config->write_file( ... );

=cut

use VRPipe::Base;

class VRPipe::Config {
    use VRPipe::Base::Configuration;
    
    has database_mysql_dbname => (
        is      => 'ro',
        isa     => 'Str',
        default => 'db_name_default_from_Config',
        documentation => 'The name of the database.'
    );
    
    has database_mysql_username => (
        is      => 'ro',
        isa     => 'Str',
        default => 'db_user_default_from_Config',
        documentation => 'The username to use when connecting to the database.'
    );
}

1;