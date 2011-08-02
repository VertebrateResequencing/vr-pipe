package TestPersistentReal;
use strict;
use warnings;
use base 'Test::DBIx::Class';
use VRPipe::Persistent::SchemaBase;

BEGIN { unshift(@INC, './modules') }
VRPipe::Persistent::SchemaBase->database_deployment('testing');

sub _initialize_schema {
    my $class = shift;
    
    my $config  = { schema_class => 'VRPipe::Persistent::Schema',
                    force_drop_table => 1,
                    keep_db => 1,
                    connect_info => [VRPipe::Persistent::SchemaBase->get_dsn, VRPipe::Persistent::SchemaBase->get_user, VRPipe::Persistent::SchemaBase->get_password] };
    
    return $class->SUPER::_initialize_schema($config);
}

1;