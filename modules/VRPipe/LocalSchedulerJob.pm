use VRPipe::Base;

class VRPipe::LocalSchedulerJob extends VRPipe::Persistent {
    has 'cmd' => (is => 'rw',
                  isa => Text,
                  traits => ['VRPipe::Persistent::Attributes']);
    
    has 'cwd' => (is => 'rw',
                  isa => Dir,
                  coerce => 1,
                  traits => ['VRPipe::Persistent::Attributes']);
    
    has 'array_size' => (is => 'rw',
                         isa => IntSQL[8],
                         traits => ['VRPipe::Persistent::Attributes'],
                         default => 1);
    
    has 'creation_time' => (is => 'rw',
                            isa => Datetime,
                            coerce => 1,
                            traits => ['VRPipe::Persistent::Attributes'],
                            default => sub { DateTime->now() });
    
    __PACKAGE__->make_persistent(has_many => [jobstates => 'VRPipe::LocalSchedulerJobState']);
}