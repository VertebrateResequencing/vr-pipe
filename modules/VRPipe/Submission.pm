use VRPipe::Base;

class VRPipe::Submission extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'job' => (is => 'rw',
                  isa => Persistent,
                  coerce => 1,
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1,
                  belongs_to => 'VRPipe::Job');
    
    has 'stepstate' => (is => 'rw',
                        isa => Persistent,
                        coerce => 1,
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_key => 1,
                        belongs_to => 'VRPipe::StepState');
    
    has 'requirements' => (is => 'rw',
                           isa => Persistent,
                           coerce => 1,
                           required => 1, # even though we're not a key
                           traits => ['VRPipe::Persistent::Attributes'],
                           belongs_to => 'VRPipe::Requirements');
    
    has 'scheduler' => (is => 'rw',
                        isa => Persistent,
                        coerce => 1,
                        required => 1,
                        builder => '_build_default_scheduler',
                        traits => ['VRPipe::Persistent::Attributes'],
                        belongs_to => 'VRPipe::Scheduler');
    
    has 'retries' => (is => 'rw',
                      isa => IntSQL[4],
                      traits => ['VRPipe::Persistent::Attributes'],
                      default => 0);
    
    has 'scheduled' => (is => 'rw',
                        isa => Datetime,
                        coerce => 1,
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_nullable => 1);
    
    has 'claim' => (is => 'rw',
                    isa => 'Bool',
                    traits => ['VRPipe::Persistent::Attributes'],
                    default => 0);
    
    has 'done' => (is => 'rw',
                   isa => 'Bool',
                   traits => ['VRPipe::Persistent::Attributes'],
                   default => 0);
    
    has 'failed' => (is => 'rw',
                     isa => 'Bool',
                     traits => ['VRPipe::Persistent::Attributes'],
                     default => 0);
    
    method _build_default_scheduler {
        return VRPipe::Scheduler->get();
    }
    
    __PACKAGE__->make_persistent();
}

1;