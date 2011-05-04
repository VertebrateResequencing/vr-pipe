use VRPipe::Base;

class VRPipe::Job extends VRPipe::Persistent {
    use DateTime;
    use Cwd;
    
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'cmd' => (is => 'rw',
                  isa => Varchar[64],
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1);
    
    has 'dir' => (is => 'rw',
                  isa => Varchar[64],
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1);
    
    has 'running' => (is => 'rw',
                      isa => 'Bool',
                      traits => ['VRPipe::Persistent::Attributes'],
                      default => 0);
    
    has 'finished' => (is => 'rw',
                       isa => 'Bool',
                       traits => ['VRPipe::Persistent::Attributes'],
                       default => 0);
    
    has 'pid' => (is => 'rw',
                  isa => IntSQL[4],
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_nullable => 1);
    
    has 'host' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_nullable => 1);
    
    has 'user' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_nullable => 1);
    
    has 'heartbeat' => (is => 'rw',
                        isa => Datetime,
                        coerce => 1,
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_nullable => 1);
    
    has 'creation_time' => (is => 'rw',
                            isa => Datetime,
                            coerce => 1,
                            traits => ['VRPipe::Persistent::Attributes'],
                            default => sub { DateTime->now() });
    
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
    
    method ok {
        if ($self->finished) {
            # check it was successful...
            return 1;
        }
        return 0;
    }
    
    __PACKAGE__->make_persistent();
    
    
    around get (ClassName|Object $self: %args) {
        unless (exists $args{dir}) {
            $args{dir} = cwd();
        }
        return $self->$orig(%args);
    }
}

1;