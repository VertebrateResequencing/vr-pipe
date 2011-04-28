use VRPipe::Base;

class VRPipe::Scheduler extends VRPipe::Persistent {
    use VRPipe::Config;
    my $vrp_config = VRPipe::Config->new();
    use VRPipe::Persistent::SchemaBase;
    
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'type' => (is => 'rw',
                   isa => Varchar[64],
                   builder => 'default_type',
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1,
                   allow_key_to_default => 1);
    
    has 'output_root' => (is => 'rw',
                          isa => Varchar[64],
                          builder => 'default_output_root',
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1,
                          allow_key_to_default => 1);
    
    method default_type (ClassName|Object $self:) {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment.'_scheduler';
        return $vrp_config->$method_name();
    }
    
    method default_output_root (ClassName|Object $self:) {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment.'_scheduler_output_root';
        return $vrp_config->$method_name();
    }
    
    __PACKAGE__->make_persistent();
}

1;