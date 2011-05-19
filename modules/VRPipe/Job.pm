use VRPipe::Base;

class VRPipe::Job extends VRPipe::Persistent {
    use DateTime;
    use Cwd;
    use Sys::Hostname;
    
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'cmd' => (is => 'rw',
                  isa => Varchar[256],
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1);
    
    has 'dir' => (is => 'rw',
                  isa => Varchar[256],
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
    
    has 'exit_code' => (is => 'rw',
                        isa => IntSQL[5],
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_nullable => 1);
    
    has 'pid' => (is => 'rw',
                  isa => IntSQL[6],
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
    
    has 'heartbeat_interval' => (is => 'rw',
                                 isa => PositiveInt,
                                 lazy => 1,
                                 builder => '_build_default_heartbeat_interval');
    
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
    
    __PACKAGE__->make_persistent();
    
    method _build_default_heartbeat_interval {
        if (VRPipe::Persistent::SchemaBase->database_deployment eq 'testing') {
            return 3;
        }
        else {
            return 60;
        }
    }
    
    around get (ClassName|Object $self: %args) {
        unless (exists $args{dir}) {
            $args{dir} = cwd();
        }
        else {
            my $dir = dir($args{dir});
            unless (-d $dir) {
                $dir->mkpath || $self->throw("job output directory '$dir' could not be created");
            }
        }
        return $self->$orig(%args);
    }
    
    method ok {
        if ($self->finished && ! $self->exit_code) {
            return 1;
        }
        return 0;
    }
    
    method pending {
        return ! ($self->running || $self->finished);
    }
    
    method run (Bool $block_and_skip_if_ok?) {
        # This sets the running state in db, then chdir to ->dir, then forks to run
        # run the ->cmd (updating ->pid, ->host, ->user and ->start_time), then when
        # finished updates ->running, ->finished and success/error-related methods as
        # appropriate. If ->run() is called when the job is not in ->pending() state,
        # it will by default return if ->finished && ->ok, otherwise will throw an
        # error. Change that behaviour:
        #$job->run(block_and_skip_if_ok => 1);
        # this waits until the job has finished running, runs it again if it failed,
        # otherwise does nothing so that $job->ok will be true. This is useful if you
        # start a bunch of tasks that all need to first index a reference file before
        # doing something on their own input files: the reference index job would be
        # run with block_and_skip_if_ok.
        # While running we update heartbeat() in the db every 30mins so that other
        # processes can query and see if we're still alive.
        
        #my $cmd = VRPipe::Cmd->new(job_id => 2445);
        #$cmd->run;
        
        # check we're allowed to run
        unless ($self->pending) {
            if ($block_and_skip_if_ok) {
                # wait until the job has finished running, run it again if it
                # failed, otherwise do nothing so that $job->ok will be true
                while (1) {
                    sleep(60);
                    if ($self->finished) {
                        if ($self->ok) {
                            return;
                        }
                        else {
                            # *** do some kind of reset on failure?
                        }
                        last;
                    }
                }
            }
            elsif ($self->ok) {
                return;
            }
            else {
                $self->throw("Job ".$self->id." could not be run because it was not in the pending state");
            }
        }
        
        #*** race condition between the pending check and claiming this to start
        #    running...
        
        # go ahead and run the cmd, also forking to run a heartbeat process
        # at the same time
        $self->running(1);
        $self->finished(0);
        $self->exit_code(undef);
        $self->update;
        
        $self->disconnect;
        my $cmd_pid = fork();
        my $heartbeat_pid;
        if (! defined $cmd_pid) {
            $self->throw("attempt to fork cmd failed: $!");
        }
        elsif ($cmd_pid) {
            # parent, start a heartbeat process to monitor the cmd
            $heartbeat_pid = fork();
            if (! defined $heartbeat_pid) {
                $self->throw("attempt to fork for heartbeat failed: $!");
            }
            elsif ($heartbeat_pid == 0) {
                # child, initiate a heartbeat that will end when the cmd stops
                # running
                my $interval = $self->heartbeat_interval;
                sleep(1);
                while (1) {
                    my $still_running = kill(0, $cmd_pid);
                    last unless $still_running;
                    $self->heartbeat(DateTime->now());
                    $self->update;
                    $self->disconnect;
                    sleep $interval;
                }
                exit(0);
            }
        }
        elsif ($cmd_pid == 0) {
            # child, actually run the command after changing to the correct dir
            # and redirecting output to files
            $self->pid($$);
            $self->host(hostname());
            $self->user(getlogin());
            $self->start_time(DateTime->now());
            $self->update;
            
            # get all info from db and disconnect before using the info below
            my $dir = $self->dir;
            my $cmd = $self->cmd;
            my $stdout_file = $self->stdout_file;
            my $stderr_file = $self->stderr_file;
            $self->disconnect;
            
            open STDOUT, '>', $stdout_file or $self->throw("Can't redirect STDOUT: $!");
            open STDERR, '>', $stderr_file or $self->throw("Can't redirect STDERR: $!");
            chdir($dir);
            exec($cmd);
        }
        
        # wait for the cmd to finish
        my $exit_code = $self->_wait_for_child($cmd_pid);
        
        # update db details for the job
        $self->end_time(DateTime->now());
        $self->running(0);
        $self->finished(1);
        $self->exit_code($exit_code);
        $self->update;
        
        # reap the heartbeat
        $self->_wait_for_child($heartbeat_pid);
    }
    
    method stdout_file {
        my $file = $self->_std_file || return;
        return file($self->dir, $file.'.o');
    }
    method stderr_file {
        my $file = $self->_std_file || return;
        return file($self->dir, $file.'.e');
    }
    method _std_file {
        my $pid = $self->pid || return;
        my $host = $self->host;
        return ".host_$host.pid_$pid";
    }
    
    method _wait_for_child (Int $pid) {
        if (waitpid($pid, 0) > 0) {
            return $?;
            #*** below code should be used elsewhere by some kind of front-end
            #    reporting method, which would access the raw exit code from
            #    ->exit_code instead of $?
            #my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
            #if ($core) {
            #    $self->warn("$identifier dumped core");
            #    return $core;
            #}
            #elsif ($sig == 9) {
            #    $self->warn("$identifier was killed");
            #    return $sig;
            #}
            #elsif ($rc) {
            #    my $message = $sig ? " after receiving signal $sig" : '';
            #    $self->warn("$identifier exited with code $rc$message");
            #    return $rc;
            #}
            #else {
            #    #warn "child $pid finished OK ($?)\n";
            #    return $rc;
            #}
        }
        else {
            my $identifier = "child process $pid for command [".$self->cmd."]";
            $self->warn("$identifier disappeared!");
        }
    }
    
    method reset {
        return if $self->running;
        #*** race condition
        $self->running(0);
        $self->finished(0);
        $self->exit_code(undef);
        $self->pid(undef);
        $self->host(undef);
        $self->user(undef);
        $self->heartbeat(undef);
        $self->creation_time(DateTime->now());
        $self->start_time(undef);
        $self->end_time(undef);
        $self->update;
        return 1;
    }
}

1;