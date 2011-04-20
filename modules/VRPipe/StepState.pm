use VRPipe::Base;

class VRPipe::StepState extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'stepmember' => (is => 'rw',
                         isa => IntSQL[16],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_key => 1);
    
    has 'dataelement' => (is => 'rw',
                          isa => IntSQL[16],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1);
    
    has 'pipelinesetup' => (is => 'rw',
                            isa => IntSQL[16],
                            traits => ['VRPipe::Persistent::Attributes'],
                            is_key => 1);
    
    has 'complete' => (is => 'rw',
                       isa => 'Bool',
                       traits => ['VRPipe::Persistent::Attributes']);
    
    __PACKAGE__->make_persistent(belongs_to => [pipelinesetup => 'VRPipe::PipelineSetup'],
                                 has_one => [dataelement => 'VRPipe::DataElement'],
                                 has_one => [stepmember => 'VRPipe::StepMember'],
                                 has_many => [submissions => 'VRPipe::Submission']);
}

1;