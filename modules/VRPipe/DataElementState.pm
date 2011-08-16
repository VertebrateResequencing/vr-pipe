use VRPipe::Base;

class VRPipe::DataElementState extends VRPipe::Persistent {
    has 'pipelinesetup' => (is => 'rw',
                            isa => Persistent,
                            coerce => 1,
                            traits => ['VRPipe::Persistent::Attributes'],
                            is_key => 1,
                            belongs_to => 'VRPipe::PipelineSetup');
    
    has 'dataelement' => (is => 'rw',
                          isa => Persistent,
                          coerce => 1,
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1,
                          belongs_to => 'VRPipe::DataElement');
    
    has 'completed_steps' => (is => 'rw',
                              isa => IntSQL[4],
                              traits => ['VRPipe::Persistent::Attributes'],
                              default => 0);
    
    __PACKAGE__->make_persistent();
}

1;