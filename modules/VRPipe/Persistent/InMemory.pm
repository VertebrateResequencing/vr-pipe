
=head1 NAME

VRPipe::Persistent::InMemory - interface to an in-memory datastore

=head1 SYNOPSIS
    
[...]

=head1 DESCRIPTION

For high-performance access to data we need right now from thousands of
processes at the same time we use an in-memory datastore that handles locking
reliably - unlike relational databases such as MySQL.

This is essentially a wrapper around Redis, providing functions needed by
VRPipe::Persistent classes.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::Persistent::InMemory {
    use VRPipe::Config;
    use VRPipe::Persistent::SchemaBase;
    use Redis;
    use Sys::Hostname;
    use Time::Format;
    use Email::Sender::Simple;
    use Email::Simple::Creator;
    use AnyEvent;
    use Scalar::Util qw(weaken isweak);
    use Time::HiRes qw(sleep);
    
    has '_maintenance_children' => (
        is      => 'ro',
        isa     => 'HashRef[ArrayRef[Int]]',
        traits  => ['Hash'],
        default => sub { {} },
        handles => {
            '_add_maintenance_child'    => 'set',
            '_have_maintenance_child'   => 'exists',
            '_delete_maintenance_child' => 'delete',
            '_all_maintenance_children' => 'values'
        },
    );
    
    our $vrp_config   = VRPipe::Config->new();
    our $email_domain = $vrp_config->email_domain();
    our $admin_email  = $vrp_config->admin_user() . '@' . $email_domain;
    our %redis_instances;
    
    our ($deployment, $reconnect_time, $log_dir, $log_file, $redis_port, $redis_server);
    
    sub BUILD {
        my $self = shift;
        
        unless (defined $redis_port) {
            $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
            $reconnect_time = $deployment eq 'production' ? 60 : 6;
            
            my $method_name = $deployment . '_logging_directory';
            $log_dir = $vrp_config->$method_name();
            $log_dir = "$log_dir";
            
            my $log_basename = VRPipe::Persistent::SchemaBase->get_dsn;
            $log_basename =~ s/\W/_/g;
            $log_file = file($log_dir, 'vrpipe-server.' . $log_basename . '.log')->stringify;
            
            $method_name = $deployment . '_redis_port';
            $redis_port  = $vrp_config->$method_name();
            $redis_port  = "$redis_port";
        }
    }
    
    sub _redis_server {
        my $self = shift;
        
        unless ($redis_server) {
            # get the hostname from our redis host file; if there isn't one
            # start the redis server first
            my $host_file = file($log_dir, 'redis.host');
            
            my $hostname;
            if (-s $host_file) {
                my $redis_host = $host_file->slurp;
                chomp($redis_host);
                my $server = $redis_host . ':' . $redis_port;
                
                # check that the redis server is still alive on that host
                my $redis;
                eval { $redis = Redis->new(server => $server, reconnect => $reconnect_time, encoding => undef); };
                if ($@) {
                    $self->log("Could not connect to Redis server \@ $server; assuming it is dead ($@)");
                    unlink($host_file);
                }
                else {
                    $hostname = $redis_host;
                    $redis->quit;
                    $self->debug("Picked up existing Redis server from host file");
                }
            }
            
            unless ($hostname) {
                # start the server
                my $redis_pid_file = file($log_dir, 'redis.pid');
                
                # check there isn't already a server running on localhost
                my $already_running = 0;
                if (-s $redis_pid_file) {
                    my $redis_pid = $redis_pid_file->slurp;
                    chomp($redis_pid);
                    if ($redis_pid) {
                        unless (kill(0, $redis_pid)) {
                            unlink($redis_pid_file);
                            $self->log("A previous Redis server had died");
                        }
                        else {
                            $self->log("A previous Redis server is still running");
                            $already_running = 1;
                        }
                    }
                }
                
                unless ($already_running) {
                    my $redis_port       = $redis_port;
                    my $redis_server_cmd = "redis-server --daemonize yes --pidfile $redis_pid_file --port $redis_port --timeout 180 --tcp-keepalive 60 --loglevel notice --logfile $log_dir/redis.log --set-max-intset-entries 2048";
                    system($redis_server_cmd);
                    sleep(1);
                    
                    unless (-s $redis_pid_file) {
                        $self->throw("Failed to start the Redis server; see $log_dir/redis.log for details");
                    }
                    else {
                        my $redis_pid = $redis_pid_file->slurp;
                        chomp($redis_pid);
                        $self->log("Started a Redis server with pid $redis_pid");
                    }
                }
                
                $hostname = hostname();
            }
            
            unless (-s $host_file) {
                my $fh = $host_file->openw;
                print $fh $hostname, "\n";
                close($fh);
                chmod 0644, $host_file;
                $self->debug("Wrote Redis host file $host_file");
            }
            
            $redis_server = $hostname . ':' . $redis_port;
        }
        
        return $redis_server;
    }
    
    method _redis {
        my $redis;
        
        # we have a class variable to store a redis instance so we can reuse
        # the connection once made, but we make sure to get a new instance for
        # each pid in case we fork - we can't use the same instance in
        # multiple processes at once
        unless (defined $redis_instances{$$}) {
            # try 3 times to get a working redis instance
            for my $try (1 .. 3) {
                eval { $redis = Redis->new(server => $self->_redis_server, reconnect => $reconnect_time, encoding => undef); };
                if ($@) {
                    if ($try == 3) {
                        $self->throw($@);
                    }
                    else {
                        sleep(1);
                    }
                }
                else {
                    last;
                }
            }
            
            $redis_instances{$$} = $redis;
        }
        else {
            $redis = $redis_instances{$$};
            # *** check it works with a ->ping? Or is that too expensive?
        }
        
        return $redis;
    }
    
    method datastore_ok {
        my $redis = $self->_redis;
        my $pong  = $redis->ping;
        return $redis->ping eq 'PONG';
    }
    
    method terminate_datastore {
        my $redis_pid_file = file($log_dir, 'redis.pid');
        if (-s $redis_pid_file) {
            my $redis_pid = $redis_pid_file->slurp;
            chomp($redis_pid);
            if ($redis_pid) {
                kill(9, $redis_pid);
                unlink($redis_pid_file);
                unlink(file($log_dir, 'redis.host'));
            }
        }
    }
    
    method lock (Str $key!, Int :$unlock_after = 900, Str :$key_prefix = 'lock', Bool :$non_exclusive = 0) {
        my $redis     = $self->_redis;
        my $redis_key = $key_prefix . '.' . $key;
        if ($non_exclusive) {
            # don't allow even a non_exclusive lock if we currently have an
            # exclusive one
            my $val = $redis->get($redis_key);
            return if $val;
        }
        
        my $got_lock = $redis->set(
            $redis_key => $non_exclusive ? 0 : hostname() . '!.!' . $$,
            EX => $unlock_after,
            $non_exclusive ? () : ('NX')
        );
        return $got_lock;
    }
    
    method _own_lock (Str $key!, Int $my_pid = $$) {
        my $redis = $self->_redis;
        my $val   = $redis->get($key);
        if ($val) {
            my ($hostname, $pid) = split('!.!', $val);
            unless ($hostname && $pid) {
                $self->debug("_own_lock failed for $key because it was set with an unexpected value of $val");
                return 0;
            }
            if ($hostname eq hostname() && $pid == $my_pid) {
                return 1;
            }
            else {
                $self->debug("_own_lock failed for $key because that is owned by $hostname|$pid, not " . hostname() . "|$my_pid");
            }
        }
        else {
            $self->debug("_own_lock failed for $key because that key is not set");
        }
        
        return 0;
    }
    
    method note (Str $key!, Int :$forget_after = 900) {
        return $self->lock($key, unlock_after => $forget_after, key_prefix => 'note', non_exclusive => 1);
    }
    
    method refresh_lock (Str $key!, Int :$unlock_after = 900, Str :$key_prefix = 'lock', :$lock_owners_pid?) {
        my $redis_key = $key_prefix . '.' . $key;
        return unless $self->_own_lock($redis_key, $lock_owners_pid || $$);
        
        my $redis = $self->_redis;
        my $refreshed = $redis->expire($redis_key, $unlock_after);
        unless ($refreshed) {
            warn "pid $$ unable to refresh redis lock for $redis_key";
            
            # presumably the redis server went down and we lost the lock;
            # if the server came back let's try and create the lock again
            $refreshed = $redis->set($redis_key => 1, EX => $unlock_after, 'NX');
        }
        
        return $refreshed;
    }
    
    method unlock (Str $key!, Str :$key_prefix = 'lock') {
        my $redis_key = $key_prefix . '.' . $key;
        return unless $self->_own_lock($redis_key);
        
        if ($self->_have_maintenance_child($redis_key)) {
            my ($owner_pid, $child_pid) = @{ $self->_delete_maintenance_child($redis_key) };
            if ($owner_pid == $$) {
                kill(9, $child_pid);
                waitpid $child_pid, 0;
            }
        }
        
        return $self->_redis->del($redis_key);
    }
    
    sub DEMOLISH {
        my $self = shift;
        foreach my $ref ($self->_all_maintenance_children) {
            my ($owner_pid, $child_pid) = @$ref;
            if ($owner_pid == $$) {
                kill(9, $child_pid);
                waitpid $child_pid, 0;
            }
        }
    }
    
    method forget_note (Str $key!) {
        my $redis_key = 'note.' . $key;
        return $self->_redis->del($redis_key);
    }
    
    method locked (Str $key!, Str :$key_prefix = 'lock', Bool :$by_me = 0) {
        my $redis_key = $key_prefix . '.' . $key;
        if ($by_me) {
            return $self->_own_lock($redis_key);
        }
        return $self->_redis->exists($redis_key);
    }
    
    method noted (Str $key!) {
        return $self->locked($key, key_prefix => 'note');
    }
    
    method block_until_unlocked (Str $key!, Int :$check_every = 2, Str :$key_prefix = 'lock') {
        while ($self->locked($key, key_prefix => $key_prefix)) {
            sleep($check_every);
        }
    }
    
    method block_until_locked (Str $key!, Int :$check_every = 2, Int :$unlock_after = 900, Str :$key_prefix = 'lock') {
        my $redis_key = $key_prefix . '.' . $key;
        return if $self->_own_lock($redis_key);
        my $sleep_time = 0.01;
        while (!$self->lock($key, unlock_after => $unlock_after, key_prefix => $key_prefix)) {
            if ($sleep_time >= $check_every) {
                $sleep_time = $check_every;
            }
            else {
                $sleep_time *= 2;
            }
            
            sleep($sleep_time);
        }
    }
    
    method maintain_lock (Str $key!, Int :$refresh_every?, Int :$leeway_multiplier?, Str :$key_prefix = 'lock') {
        my $redis_key = $key_prefix . '.' . $key;
        $self->throw("maintain_lock() cannot be used unless we own the lock") unless $self->_own_lock($redis_key);
        
        unless ($self->_have_maintenance_child($redis_key)) {
            unless ($refresh_every) {
                $refresh_every = $deployment eq 'testing' ? 3 : 60;
            }
            
            unless ($leeway_multiplier) {
                if ($refresh_every <= 60) {
                    $leeway_multiplier = 15;
                }
                elsif ($refresh_every <= 300) {
                    $leeway_multiplier = 3;
                }
                else {
                    $leeway_multiplier = 2;
                }
            }
            my $survival_time = $refresh_every * $leeway_multiplier;
            
            # we used to just have an AnyEvent->timer that called refresh_lock
            # periodically, but that timer doesn't fire while normal code is
            # busy running. I also couldn't get a 'cede' type method working
            # due to recursion in AnyEvent watchers. Instead we fork to
            # periodically calls refresh_lock
            # *** should we use Proc::FastSpawn or AnyEvent::Fork for greater
            # speed and lower memory?
            my $my_pid   = $$;
            my $lock_pid = fork();
            if (!defined $lock_pid) {
                $self->throw("attempt to fork for lock failed: $!");
            }
            elsif ($lock_pid == 0) {
                # child, initiate a lock that will end when the parent stops
                # running
                while (1) {
                    kill(0, $my_pid) || last;
                    last unless $self->_own_lock($redis_key, $my_pid);
                    $self->refresh_lock($key, key_prefix => $key_prefix, unlock_after => $survival_time, lock_owners_pid => $my_pid);
                    sleep $refresh_every;
                }
                exit(0);
            }
            
            $self->_add_maintenance_child($redis_key, [$$, $lock_pid]);
        }
        
        return 1;
    }
    
    method enqueue (Str $key!, Str $value!) {
        return $self->_redis->sadd('queue.' . $key, $value);
    }
    
    method dequeue (Str $key!, ArrayRef $members?) {
        my $redis_key = 'queue.' . $key;
        if ($members) {
            return $self->_redis->srem($redis_key, @$members);
        }
        else {
            return $self->_redis->spop($redis_key);
        }
    }
    
    method queue (Str $key!) {
        return $self->_redis->smembers('queue.' . $key);
    }
    
    method drop_queue (Str $key!) {
        return $self->_redis->del('queue.' . $key);
    }
    
    method log_stderr {
        my $ok;
        eval {
            # *** for some unknown reason, if I use $self->_redis() here the
            # server breaks due to "corrupt" redis instances... I don't get it!!
            my $redis = Redis->new(server => $self->_redis_server, reconnect => $reconnect_time, encoding => undef);
            $ok = $redis->ping;
        };
        
        if ($ok) {
            tie *STDERR, 'VRPipe::Persistent::InMemory', $self->_redis_server, $reconnect_time, $log_file;
        }
        else {
            $ok = open(STDERR, '>>', $log_file);
        }
    }
    
    sub TIEHANDLE {
        my ($pkg, $server, $reconnect, $log_file) = @_;
        return bless { server => $server, reconnect => $reconnect, log_file => $log_file }, $pkg;
    }
    
    sub PRINT {
        my $config = shift;
        my ($ok, $redis);
        eval {
            $redis = Redis->new(server => $config->{server}, reconnect => $config->{reconnect}, encoding => undef);
            $ok = $redis->ping;
        };
        
        if ($ok) {
            $redis->rpush('stderr', @_);
        }
        else {
            open(my $fh, '>>', $config->{log_file}) || return;
            print $fh @_;
            close($fh);
        }
    }
    
    method write_stderr {
        open(my $fh, '>>', $log_file) || return;
        my $redis = $self->_redis || return;
        while (my $line = $redis->lpop('stderr')) {
            print $fh $line;
        }
        close($fh);
    }
    
    method log (Str $msg!, ArrayRef[Str] :$email_to?, Bool :$email_admin?, Str :$subject?, Str :$long_msg?, Bool :$force_when_testing?) {
        chomp($msg);
        
        # we'll just warn, which will end up in the log file automatically; we
        # do not try and do anything fancy with flock() since that can be too
        # slow or end up locking the log file forever
        my $log_msg = "$time{'yyyy-mm-dd hh:mm:ss'} | pid $$ | $msg\n";
        warn $log_msg;
        
        if (($force_when_testing || $deployment eq 'production') && ($email_to || $email_admin)) {
            # email the desired users
            my $email = Email::Simple->create(
                header => [
                    To => $email_to ? join(', ', map { "$_\@$email_domain" } @$email_to) : $admin_email,
                    $email_admin && $email_to ? (Cc => $admin_email) : (),
                    From    => qq["VRPipe Server" <$admin_email>],
                    Subject => $subject || "VRPipe Server message",
                ],
                body => $msg . "\n" . ($long_msg || ''),
            );
            my $sent = Email::Sender::Simple->try_to_send($email);
            
            unless ($sent) {
                warn "$time{'yyyy-mm-dd hh:mm:ss'}: previous message failed to get sent to [", join(', ', ($email_to ? @$email_to : (), $email_admin ? '' : ())), "]\n";
            }
        }
    }
    
    method debug (Str $msg!) {
        return unless $self->verbose > 0;
        return $self->log($msg);
    }
}

1;
