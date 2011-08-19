use VRPipe::Base;

class VRPipe::Submission extends VRPipe::Persistent {
    use Devel::GlobalDestruction;
    use DateTime;
    use VRPipe::Parser;
    
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
        $self->throw("Cannot call update_status when the job ".$self->job->id." is not finished") unless $self->job->finished;
        return if $self->done || $self->failed;
        
        if ($self->job->ok) {
            $self->_done(1);
            $self->_failed(0);
        }
        else {
            $self->_done(0);
            $self->_failed(1);
        }
        
        unless ($self->_hid) {
            # we're a submission for a job that completed in a different
            # submission, so we didn't actually run and have no output
            $self->update;
            return;
        }
        
        #*** these 2 calls are probably the cause of massive delays when going
        #    from jobs being finished to submissions being done... can we
        #    optimise?
        #    It is actualy the archive_output that causes the ~1s delays,
        #    and it is not the filesystem operations, but the retrieval of
        #    the 8 VRPipe::File objects... at the very least path column of file
        #    table must be indexed, or delays go up to ~6s+
        my $t1 = time();
        $self->sync_scheduler;
        my $t2 = time();
        $self->archive_output;
        my $t3 = time();
        $self->_sid(undef);
        if ($t3 - $t1 > 2) {
            my $ss_t = $t2 - $t1;
            my $ao_t = $t3 - $t2;
            $self->debug("in update_status, sync_scheduler took $ss_t seconds; archive_output took $ao_t seconds");
        }
        
        $self->update;
    }
    
    method sync_scheduler {
        my $sid = $self->sid;
        unless ($sid) {
            return;
        }
        
        # assume that if its been more than a minute since the job ended, the
        # scheduler must have generated its output by now, so we don't have to
        # ask the scheduler about this sid
        my $job = $self->job;
        my $end_time = $self->job->end_time;
        if ($end_time && (time() - $end_time->epoch > 60)) {
            return;
        }
        
        $self->scheduler->wait_for_sid($sid, $self->_aid, 5);
    }
    
    method archive_output {
        my $jso = $self->job->stdout_file || $self->warn("no job stdout_file for job ".$self->job->id);
        return unless $jso;
        my $sso = $self->job_stdout_file;
        $self->move($jso, $sso) if $jso->e;
        my $jse = $self->job->stderr_file;
        my $sse = $self->job_stderr_file;
        $self->move($jse, $sse) if $jse->e;
        
        my $scso = $self->scheduler_stdout_file(orig => 1);
        my $dest = $self->scheduler_stdout_file;
        $self->move($scso, $dest) if $scso->e;
        my $scse = $self->scheduler_stderr_file(orig => 1);
        $dest = $self->scheduler_stderr_file;
        $self->move($scse, $dest) if $scse->e;
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
    
    sub DEMOLISH {
        return if in_global_destruction;
        my $self = shift;
        $self->release if $self->in_storage;
    }
    
    __PACKAGE__->make_persistent();
    
    # where have our scheduler and job stdout/err files gone?
    method std_dir {
        my ($for) = $self->_for;
        $for || return;
        return $self->scheduler->output_dir($for);
    }
    
    method _for {
        my $hid = $self->_hid || return;
        my ($for, $index);
        if ($self->_aid) {
            $for = VRPipe::PersistentArray->get(id => $hid);
            $index = $self->_aid;
        }
        else {
            $for = $self;
        }
        return ($for, $index);
    }
    
    method _scheduler_std_file (Str $method where {$_ eq 'scheduler_output_file' || $_ eq 'scheduler_error_file'}, Str $type where {$_ eq 'lsf' || $_ eq 'txt'}, Bool $orig = 0) {
        my $std_dir = $self->std_dir || return;
        my $std_io_file = $self->scheduler->$method($std_dir);
        my (undef, $index) = $self->_for;
        if ($index) {
            $std_io_file .= '.'.$index;
        }
        $std_io_file .= '.r'.$self->retries unless $orig;
        return VRPipe::File->get(path => $std_io_file, type => $type);
    }
    method scheduler_stdout_file (Bool :$orig = 0) {
        return $self->_scheduler_std_file('scheduler_output_file', 'lsf', $orig); #*** should not be hard-coded for lsf
    }
    method scheduler_stderr_file (Bool :$orig = 0) {
        return $self->_scheduler_std_file('scheduler_error_file', 'txt', $orig);
    }
    method scheduler_stdout {
        my $file = $self->scheduler_stdout_file;
        $file->s || return;
        return VRPipe::Parser->create('lsf', {file => $file}); #*** should not be hard-coded for lsf
    }
    method scheduler_stderr {
        my $file = $self->scheduler_stderr_file;
        $file->s || return;
        return $file->slurp;
    }
    
    method _job_std_file (Str $kind where {$_ eq 'out' || $_ eq 'err'}) {
        # this is where we want the job stdo/e to be archived to, not where the
        # job initially spits it out to
        my $std_dir = $self->std_dir || return;
        my (undef, $index) = $self->_for;
        if ($index) {
            $index = '.'.$index;
        }
        else {
            $index = '';
        }
        return VRPipe::File->get(path => file($std_dir, 'job_std'.$kind.$index.'.r'.$self->retries), type => 'txt');
    }
    method job_stdout_file {
        return $self->_job_std_file('out');
    }
    method job_stderr_file {
        return $self->_job_std_file('err');
    }
    method job_stdout {
        my $file = $self->job_stdout_file;
        $file->s || return;
        return $file->slurp;
    }
    method job_stderr {
        my $file = $self->job_stderr_file;
        $file->s || return;
        return $file->slurp;
    }
    
    method pend_time {
        my $job = $self->job;
        my $scheduled_time = $self->_scheduled || $self->throw("called pend_time, yet submission ".$self->id." has not been scheduled!");
        $scheduled_time = $scheduled_time->epoch;
        if ($job->running || $job->finished) {
            return $job->start_time->epoch -$scheduled_time;
        }
        else {
            return time() - $scheduled_time;
        }
    }
    
    method unschedule_if_not_pending {
        my $sid = $self->sid || return;
        my $aid = $self->_aid;
        my $scheduler = $self->scheduler;
        
        my $status = $scheduler->sid_status($sid, $aid);
        return if $status eq 'PEND';
        $scheduler->kill_sid($sid, $aid, 5);
        
        $self->_reset;
    }
    
    method retry {
        return unless ($self->done || $self->failed);
        
        # reset the job
        $self->_reset_job;
        
        # reset ourself
        my $retries = $self->retries;
        $self->retries($retries + 1);
        $self->_reset;
    }
    
    method start_over {
        # reset the job
        $self->_reset_job;
        
        # reset outself and also set retries to 0
        $self->retries(0);
        $self->_reset;
    }
    
    method _reset_job {
        my $job = $self->job;
        unless ($job->finished) {
            $job->kill_job;
        }
        
        my $ofiles = $job->output_files;
        foreach my $file (@$ofiles) {
            $file->unlink;
        }
        
        $job->reset_job;
    }
    
    method _reset {
        $self->_sid(undef);
        $self->_failed(0);
        $self->_done(0);
        $self->_aid(undef);
        $self->_scheduled(undef);
        $self->_claim(0);
        $self->_hid(undef);
        $self->update;
    }
}

1;