
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

Copyright (c) 2013-2015 Genome Research Limited.

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
    use Bytes::Random::Secure;
    use VRPipe::Interface::BackEnd;
    
    our $vrp_config   = VRPipe::Config->new();
    our $email_domain = $vrp_config->email_domain();
    our $admin_email  = $vrp_config->admin_user() . '@' . $email_domain;
    our %redis_instances;
    
    our ($deployment, $reconnect_time, $log_dir, $log_file, $redis_port, $redis_server, $backend);
    
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
            
            $backend = VRPipe::Interface::BackEnd->new(deployment => $deployment);
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
    
    # default of unlock_after => 0 means the lock is held "forever", but any
    # other process trying to get a lock will see if this process is still
    # alive, and if not, the initial lock will be lost
    method lock (Str $key!, Int :$unlock_after = 0, Str :$key_prefix = 'lock', Bool :$non_exclusive = 0, Bool :$debug = 0, Bool :$check_via_ssh = 1) {
        my $redis     = $self->_redis;
        my $redis_key = $key_prefix . '.' . $key;
        
        my $valid_lock = $self->_check_existing_lock($redis_key, debug => $debug, check_via_ssh => $check_via_ssh);
        
        if ($non_exclusive) {
            # don't allow even a non_exclusive lock if we currently have an
            # exclusive one
            return if $valid_lock;
        }
        
        my $redis_key_value = $non_exclusive ? 0 : hostname() . '!.!' . $$;
        
        my $got_lock = $redis->set(
            $redis_key => $redis_key_value,
            $unlock_after ? (EX => $unlock_after) : (),
            $non_exclusive ? () : ('NX')
        );
        
        if ($debug) {
            if ($got_lock) {
                warn " lock for $redis_key => $redis_key_value worked\n";
            }
            else {
                warn " lock for $redis_key => $redis_key_value failed because $redis_key is currently set to $valid_lock\n";
            }
        }
        
        return $got_lock;
    }
    
    method _own_lock (Str $key!, Int $my_pid = $$) {
        my $redis = $self->_redis;
        my $val   = $self->_check_existing_lock($key);
        if (defined $val) {
            my ($hostname, $pid) = split('!.!', $val);
            if ($hostname eq hostname() && $pid == $my_pid) {
                return 1;
            }
            else {
                $self->debug("_own_lock failed for $key because that is owned by $hostname|$pid, not " . hostname() . "|$my_pid");
            }
        }
        return 0;
    }
    
    method _check_existing_lock (Str $key!, Bool :$debug = 0, Bool :$check_via_ssh = 1) {
        my $redis = $self->_redis;
        my $val   = $redis->get($key);
        
        if ($val) {
            my ($hostname, $pid) = split('!.!', $val);
            unless ($hostname && $pid) {
                warn " lock for $key had an invalid value of $val, so removing it\n" if $debug;
                $redis->del($key);
                return;
            }
            
            # check the current lock is for a currently living pid, or has a
            # specified EX
            my $ttl = $redis->ttl($key);
            unless ($ttl && $ttl > 0) {
                if ($hostname eq hostname()) {
                    unless (kill(0, $pid)) {
                        warn " lock for $key was initiated on host $hostname but the pid $pid is dead, so removing it\n" if $debug;
                        $redis->del($key);
                        return;
                    }
                }
                elsif ($check_via_ssh) {
                    my $return;
                    eval {
                        local $SIG{ALRM} = sub { die "alarm\n" };
                        alarm 15;
                        # we use perl for testing the pid because don't know
                        # what shell the user has and what syntax to use
                        $return = $backend->ssh($hostname, "perl -e 'print qq[pid_alive] if kill(0, $pid)'");
                        alarm 0;
                    };
                    
                    if (!$return || $return !~ /pid_alive/) {
                        warn " lock for $key was initiated on host $hostname but the pid $pid is dead, so removing it\n" if $debug;
                        $redis->del($key);
                        return;
                    }
                }
            }
        }
        
        return $val;
    }
    
    method note (Str $key!, Int :$forget_after = 900) {
        return $self->lock($key, unlock_after => $forget_after, key_prefix => 'note', non_exclusive => 1);
    }
    
    method unlock (Str $key!, Str :$key_prefix = 'lock') {
        my $redis_key = $key_prefix . '.' . $key;
        return unless $self->_own_lock($redis_key);
        return $self->_redis->del($redis_key);
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
        #else
        my $val = $self->_check_existing_lock($redis_key);
        return defined $val ? 1 : 0;
    }
    
    method noted (Str $key!) {
        return $self->locked($key, key_prefix => 'note');
    }
    
    method block_until_unlocked (Str $key!, Int :$check_every = 2, Str :$key_prefix = 'lock') {
        while ($self->locked($key, key_prefix => $key_prefix)) {
            sleep($check_every);
        }
    }
    
    method block_until_locked (Str $key!, Int :$check_every = 2, Str :$key_prefix = 'lock', Bool :$debug = 0, Int :$unlock_after = 0) {
        my $redis_key = $key_prefix . '.' . $key;
        return if $self->_own_lock($redis_key);
        my $sleep_time = 0.01;
        my $t          = time();
        while (!$self->lock($key, unlock_after => $unlock_after, key_prefix => $key_prefix, debug => $debug, check_via_ssh => time() == $t)) {
            if ($sleep_time >= $check_every) {
                $sleep_time = $check_every;
            }
            else {
                $sleep_time *= 2;
            }
            
            warn " - didn't get a lock on $redis_key, will try again in $sleep_time seconds\n" if $debug;
            sleep($sleep_time);
            
            # we'll check_via_ssh only ~once per minute
            if (time() > $t + 60) {
                $t = time();
            }
        }
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
    
    method create_session (HashRef $data!, Int :$idle_expiry = 86400, Int :$max_life = 432000) {
        # (the idle_expiry and max_life are very lax by default, designed to
        #  allow someone to keep a session for a working week if they use the
        #  session at least once per day during work hours)
        my $redis     = $self->_redis;
        my $random    = Bytes::Random::Secure->new(Bits => 128, NonBlocking => 1);
        my $key       = $random->bytes_hex(16);
        my $redis_key = 'session.' . $key;
        if ($redis->exists($redis_key)) {
            # wow, this should be basically impossible
            $self->throw("A session with key '$key' already exists, we've lost randomness!");
        }
        
        my $absolute_expiry = time() + $max_life;
        my $set             = $redis->hmset(
            $redis_key,
            %$data,
            'session.idle_expiry'     => $idle_expiry,
            'session.absolute_expiry' => $absolute_expiry
        );
        
        if ($set) {
            $redis->expire($redis_key, $idle_expiry);
            # (exireat overrides expire, so I implement my own in get_session())
            return $key;
        }
        return;
    }
    
    method get_session (Str $key!) {
        my $redis     = $self->_redis;
        my $redis_key = 'session.' . $key;
        
        my %hash = $redis->hgetall($redis_key);
        if (keys %hash) {
            # refresh expiry
            my $absolute_expiry = delete $hash{'session.absolute_expiry'};
            if (time() > $absolute_expiry) {
                $self->drop_session($key);
                return;
            }
            my $idle_expiry = delete $hash{'session.idle_expiry'};
            $redis->expire($redis_key, $idle_expiry);
            
            return \%hash;
        }
        return;
    }
    
    method session_set (Str $key!, Str $set_key!, Str $value!) {
        my $redis     = $self->_redis;
        my $redis_key = 'session.' . $key;
        if ($redis->exists($redis_key)) {
            # we only set new key vals for existing sessions
            return $redis->hset($redis_key, $set_key => $value);
        }
        return 0;
    }
    
    method session_get (Str $key!, Str $get_key!) {
        return $self->_redis->hget('session.' . $key, $get_key);
    }
    
    method session_del (Str $key!, Str $del_key!) {
        return $self->_redis->hdel('session.' . $key, $del_key);
    }
    
    method drop_session (Str $key!) {
        return $self->_redis->del('session.' . $key);
    }
    
    method rate_limit (Str $key!, Int :$per_second = 1, Bool :$punish_excess = 0) {
        my $redis     = $self->_redis;
        my $redis_key = 'rate_limit.' . $key;
        my $count     = $redis->incr($redis_key);
        
        if ($count) {
            my $expire_time = 1;
            if ($punish_excess) {
                # instead of rate limiting to 1 per second, if rate_limit is
                # called before expiration, we ramp up the expiration time
                if ($count <= 10) {
                    $expire_time = $count;
                }
                else {
                    $expire_time = 45;
                }
            }
            $redis->expire($redis_key, $expire_time);
            return $count <= $per_second ? 1 : 0;
        }
        else {
            # something went wrong in redis, treat it as limit
            return 0;
        }
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
