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
                          isa => IntSQL[16],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1);
    
    has 'pipeline' => (is => 'rw',
                       isa => IntSQL[16],
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1);
    
    has 'options' => (is => 'rw',
                      isa => Varchar[64],
                      traits => ['VRPipe::Persistent::Attributes'],
                      is_key => 1);
    
    has 'description' => (is => 'rw',
                         isa => Varchar[64],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    __PACKAGE__->make_persistent(has_one => [datasource => 'VRPipe::DataSource'],
                                 has_one => [pipeline => 'VRPipe::Pipeline'],
                                 has_many => [states => 'VRPipe::StepState']);
}

1;