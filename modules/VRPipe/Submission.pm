use VRPipe::Base;

class VRPipe::Submission extends VRPipe::Persistent {
    use DateTime;
    
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
                           # handles => [qw(memory time cpu tmp_space local_space custom)]), *** doesn't work for some reason, and we need them read-only anyway
                           belongs_to => 'VRPipe::Requirements');
    
    has 'scheduler' => (is => 'rw',
                        isa => Persistent,
                        coerce => 1,
                        required => 1,
                        builder => '_build_default_scheduler',
                        traits => ['VRPipe::Persistent::Attributes'],
                        belongs_to => 'VRPipe::Scheduler');
    
    has '_sid' => (is => 'rw',
                   isa => IntSQL[8],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_nullable => 1);
    
    has '_hid' => (is => 'rw',
                   isa => IntSQL[8],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_nullable => 1);
    
    has '_aid' => (is => 'rw',
                   isa => IntSQL[8],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_nullable => 1);
    
    has 'retries' => (is => 'rw',
                      isa => IntSQL[4],
                      traits => ['VRPipe::Persistent::Attributes'],
                      default => 0);
    
    has '_scheduled' => (is => 'rw',
                         isa => Datetime,
                         coerce => 1,
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    has '_claim' => (is => 'rw',
                     isa => 'Bool',
                     traits => ['VRPipe::Persistent::Attributes'],
                     default => 0);
    
    has '_own_claim' => (is => 'rw',
                         isa => 'Bool',
                         default => 0);
    
    has '_done' => (is => 'rw',
                    isa => 'Bool',
                    traits => ['VRPipe::Persistent::Attributes'],
                    default => 0);
    
    has '_failed' => (is => 'rw',
                      isa => 'Bool',
                      traits => ['VRPipe::Persistent::Attributes'],
                      default => 0);
    
    method _build_default_scheduler {
        return VRPipe::Scheduler->get();
    }
    
    # public getters for our private attributes
    method sid (PositiveInt $sid?) {
        if ($sid) {
            return unless $self->claim;
            
            $self->_sid($sid);
            
            $self->_scheduled(DateTime->now);
            $self->release;
            $self->update;
            
            return $sid;
        }
        else {
            return $self->_sid;
        }
    }
    
    method scheduled {
        return $self->_scheduled;
    }
    
    method done {
        return $self->_done;
    }
    
    method failed {
        return $self->_failed;
    }
    
    method claim {
        return 0 if $self->scheduled;
        return 0 if $self->sid;
        
        if ($self->_claim) {
            return $self->_own_claim ? 1 : 0;
        }
        else {
            $self->_claim(1);
            $self->update;
            $self->_own_claim(1);
            return 1;
        }
    }
    
    # scheduling-related behaviour
    method release {
        $self->_claim(0);
        $self->_own_claim(0);
        $self->update;
    }
    
    method submit {
        $self->scheduler->submit(submission => $self);
    }
    
    around scheduler {
        if ($self->scheduled) {
            return $self->$orig;
        }
        else {
            return $self->$orig(@_);
        }
    }
    
    method update_status {
        $self->throw("Cannot call update_status when the job is not finished") unless $self->job->finished;
        
        if ($self->job->ok) {
            $self->_done(1);
            $self->_failed(0);
        }
        else {
            $self->_done(0);
            $self->_failed(1);
        }
        
        #$self->sync_scheduler;
        #$self->archive_output; #*** cat the raw stdout and stderr files created by Job named after the pid in the Job working dir on to the end of files in hashed global output_dir directory
        $self->_sid(undef);
        
        $self->update;
    }
    
    # requirement passthroughs and extra_* methods
    method _add_extra (Str $type, Int $extra) {
        my $new_req = $self->requirements->clone($type => $self->$type() + $extra);
        $self->requirements($new_req);
        $self->update;
    }
    
    method memory {
        return $self->requirements->memory;
    }
    method extra_memory (Int $extra = 1000) {
        $self->_add_extra('memory', $extra);
    }
    
    method time {
        return $self->requirements->time;
    }
    method extra_time (Int $extra = 3600) {
        $self->_add_extra('time', $extra);
    }
    
    method cpu {
        return $self->requirements->cpu;
    }
    method extra_cpu (Int $extra = 1) {
        $self->_add_extra('cpu', $extra);
    }
    
    method tmp_space {
        return $self->requirements->tmp_space;
    }
    method extra_tmp_space (Int $extra = 4000) {
        $self->_add_extra('tmp_space', $extra);
    }
    
    method local_space {
        return $self->requirements->local_space;
    }
    method extra_local_space (Int $extra = 4000) {
        $self->_add_extra('local_space', $extra);
    }
    
    method custom {
        return $self->requirements->custom;
    }
    
    method DEMOLISH {
        $self->release if $self->in_storage;
    }
    
    __PACKAGE__->make_persistent();
}

1;