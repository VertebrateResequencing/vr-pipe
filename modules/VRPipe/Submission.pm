use VRPipe::Base;

class VRPipe::Submission extends VRPipe::Persistent {
    use VRPipe::Config;
    my $vrp_config = VRPipe::Config->new();
    
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'job' => (is => 'rw',
                  isa => Persistent,
                  coerce => 1,
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1);
    
    has 'stepstate' => (is => 'rw',
                        isa => IntSQL[16],
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_key => 1);
    
    has 'requirments' => (is => 'rw',
                          isa => Persistent,
                          coerce => 1,
                          required => 1,
                          traits => ['VRPipe::Persistent::Attributes']);
    
    has 'scheduler' => (is => 'rw',
                        isa => Persistent,
                        coerce => 1,
                        required => 1,
                        builder => '_build_default_scheduler',
                        traits => ['VRPipe::Persistent::Attributes']);
    
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
        my $method_name = $VRPipe::Persistent::SchemaBase::DATABASE_DEPLOYMENT.'_scheduler';
        my $scheduler = $vrp_config->$method_name();
        return VRPipe::Scheduler->get(module => "VRPipe::Schedulers::$scheduler");
    }
    
    __PACKAGE__->make_persistent(belongs_to => [stepstate => 'VRPipe::StepState'],
                                 has_one => [job => 'VRPipe::Job'],
                                 has_one => [requirements => 'VRPipe::Requirements'],
                                 has_one => [scheduler => 'VRPipe::Scheduler']);
}

1;