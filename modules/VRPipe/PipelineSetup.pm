use VRPipe::Base;

class VRPipe::PipelineSetup extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'datasource' => (is => 'rw',
                         isa => Persistent,
                         coerce => 1,
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_key => 1,
                         belongs_to => 'VRPipe::DataSource');
    
    has 'pipeline' => (is => 'rw',
                       isa => Persistent,
                       coerce => 1,
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       belongs_to => 'VRPipe::Pipeline');
    
    has 'options' => (is => 'rw',
                      isa => Varchar[64],
                      traits => ['VRPipe::Persistent::Attributes'],
                      default => '',
                      allow_key_to_default => 1,
                      is_key => 1);
    
    has 'description' => (is => 'rw',
                         isa => Varchar[64],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    __PACKAGE__->make_persistent(has_many => [states => 'VRPipe::StepState']);
}

1;