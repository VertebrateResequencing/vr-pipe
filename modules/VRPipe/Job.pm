use VRPipe::Base;

class VRPipe::Job extends VRPipe::Persistent {
    use DateTime;
    use Cwd;
    use Sys::Hostname;
    use Net::SSH qw(ssh);
    
    has 'cmd' => (is => 'rw',
                  isa => Text,
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1);
    
    has 'dir' => (is => 'rw',
                  isa => Dir,
                  coerce => 1,
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1);
    
    has 'block_and_skip_if_ok' => (is => 'rw',
                                   isa => 'Bool',
                                   traits => ['VRPipe::Persistent::Attributes'],
                                   default => 0);
    
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
    
    has 'output_files' => (is => 'rw',
                           isa => ArrayRefOfPersistent,
                           traits => ['VRPipe::Persistent::Attributes'],
                           default => sub { [] });
    
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
    
    method run (VRPipe::StepState :$stepstate?) {
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
        # While running we update heartbeat() in the db every min so that other
        # processes can query and see if we're still alive.
        
        # check we're allowed to run, in a transaction to avoid race condition
        my $schema = $self->result_source->schema;
        $schema->txn_do(sub {
            unless ($self->pending) {
                if ($self->block_and_skip_if_ok) {
                    # wait until the job has finished running, run it again if it
                    # failed, otherwise do nothing so that $job->ok will be true
                    while (1) {
                        $self->disconnect;
                        sleep(60);
                        if ($self->finished) {
                            if ($self->ok) {
                                return;
                            }
                            else {
                                # *** do some kind of reset on failure?
                                $self->throw("blocking and skipping if ok, but finished and failed... don't know what to do!");
                            }
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
            
            # go ahead and run the cmd, also forking to run a heartbeat process
            # at the same time
            $self->running(1);
            $self->start_time(DateTime->now());
            $self->finished(0);
            $self->exit_code(undef);
            $self->update;
        });
        
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
            $self->user(getlogin || getpwuid($<));
            $self->update;
            
            # get all info from db and disconnect before using the info below
            my $dir = $self->dir;
            my $cmd = $self->cmd;
            my $stdout_file = $self->stdout_file->path; #*** call here in child, then below in parent, creating 2 identical rows in db?... should not be possible!
            $stdout_file->dir->mkpath;
            my $stderr_file = $self->stderr_file->path; #*** but only 1 copy of the e?!
            $self->disconnect;
            
            open STDOUT, '>', $stdout_file or $self->throw("Can't redirect STDOUT to '$stdout_file': $!");
            open STDERR, '>', $stderr_file or $self->throw("Can't redirect STDERR to '$stderr_file': $!");
            chdir($dir);
            exec($cmd);
        }
        
        # wait for the cmd to finish
        my $exit_code = $self->_wait_for_child($cmd_pid);
        
        # update db details for the job
        $self->end_time(DateTime->now());
        $self->stdout_file->update_stats_from_disc(retries => 3);
        $self->stderr_file->update_stats_from_disc(retries => 3);
        my $o_files = $self->output_files;
        if (@$o_files) {
            foreach my $o_file (@$o_files) {
                $o_file->update_stats_from_disc(retries => 3);
            }
        }
        elsif ($stepstate) {
            $stepstate->update_output_file_stats;
        }
        $self->running(0);
        $self->finished(1);
        $self->exit_code($exit_code);
        $self->update;
        
        kill(9, $heartbeat_pid);
        # reap the heartbeat
        $self->_wait_for_child($heartbeat_pid);
    }
    
    method stdout_file {
        my $dir = $self->dir;
        my $file = $self->_std_file || return;
        return VRPipe::File->get(path => file($dir, $file.'.o'), type => 'txt');
    }
    method stderr_file {
        my $dir = $self->dir;
        my $file = $self->_std_file || return;
        return VRPipe::File->get(path => file($dir, $file.'.e'), type => 'txt');
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
    
    method time_since_heartbeat {
        return unless $self->running;
        my $heartbeat = $self->heartbeat || $self->start_time || $self->throw("job ".$self->id." is running, yet has neither a heartbeat nor a start time");
        my $t = time();
        return $t - $heartbeat->epoch;
    }
    
    method unresponsive {
        my $interval = $self->heartbeat_interval;
        my $elapsed = $self->time_since_heartbeat || 0;
        
        # try and allow for mysql server time being different to our host time,
        # and other related timing vagueries
        my $multiplier;
        if ($interval < 60) {
            $multiplier = 10;
        }
        elsif ($interval < 300) {
            $multiplier = 3;
        }
        else {
            $multiplier = 1;
        }
        
        $self->debug("   -- checking if unresponsive with a heartbeat interval of $interval, elapsed time of $elapsed s since last heartbeat, and multiplier $multiplier");
        
        return $elapsed > ($interval * $multiplier) ? 1 : 0;
    }
    
    method kill_job {
        return unless $self->running;
        my ($user, $host, $pid) = ($self->user, $self->host, $self->pid);
        if ($user && $host && $pid) {
            $self->disconnect;
            ssh("$user\@$host", "kill -9 $pid"); #*** we will fail to login with key authentication if user has never logged into this host before, and it asks a question...
                                                 #    Net::SSH::Perl is able to always log us in, but can take over a minute!
            # *** do we care if the kill fails?...
        }
        $self->running(0);
        $self->finished(1);
        $self->end_time(DateTime->now());
        $self->exit_code(9);
        $self->update;
    }
    
    method reset_job {
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