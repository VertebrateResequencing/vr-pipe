use VRPipe::Base;

class VRPipe::LocalSchedulerJobState extends VRPipe::Persistent {
    has 'localschedulerjob' => (is => 'rw',
                                isa => IntSQL[8],
                                coerce => 1,
                                belongs_to => 'VRPipe::LocalSchedulerJob',
                                traits => ['VRPipe::Persistent::Attributes'],
                                is_key => 1);
    
    has 'aid' => (is => 'rw',
                  isa => IntSQL[8],
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1);
    
    has 'start_time' => (is => 'rw',
                         isa => Datetime,
                         coerce => 1,
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    has 'end_time' => (is => 'rw',
                       isa => Datetime,
                       coerce => 1,
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_nullable => 1);
    
    has 'exit_code' => (is => 'rw',
                        isa => IntSQL[5],
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_nullable => 1);
    
    __PACKAGE__->make_persistent();
}