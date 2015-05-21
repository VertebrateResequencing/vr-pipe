
=head1 NAME

VRPipe::Interface::BackEnd - methods for the server to handle user interaction

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2014 Genome Research Limited.

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

class VRPipe::Interface::BackEnd {
    use Perl6::Form;
    use Module::Find;
    use VRPipe::Persistent::SchemaBase;
    use VRPipe::Config;
    use AnyEvent::Util qw(fork_call);
    use Twiggy::Server::TLS;
    use Plack::Request;
    use POSIX qw(setsid setuid);
    use Fcntl ':mode';
    use Cwd qw(chdir getcwd);
    use Module::Find;
    use VRPipe::Persistent::InMemory;
    
    my %extension_to_content_type = (
        gif  => 'image/gif',
        jpg  => 'image/jpeg',
        png  => 'image/png',
        csv  => 'text/csv',
        html => 'text/html',
        js   => 'text/javascript',
        css  => 'text/css',
        txt  => 'text/plain',
        zip  => 'application/zip',
        gzip => 'application/gzip',
        pdf  => 'application/pdf',
        bin  => 'application/octet-stream'
    );
    
    has 'deployment' => (
        is       => 'ro',
        isa      => 'Str',
        required => 1
    );
    
    has 'farm' => (
        is  => 'ro',
        isa => 'Str'
    );
    
    has 'port' => (
        is     => 'ro',
        isa    => PositiveInt,
        writer => '_set_port'
    );
    
    has 'tls_key' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_tls_key'
    );
    
    has 'tls_cert' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_tls_cert'
    );
    
    has 'dsn' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_dsn'
    );
    
    has 'neo4j_server_url' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_neo4j_server_url'
    );
    
    has 'scheduler' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_scheduler'
    );
    
    has 'schema' => (
        is      => 'ro',
        isa     => 'VRPipe::Persistent::Schema',
        lazy    => 1,
        builder => '_build_schema'
    );
    
    has 'umask' => (
        is     => 'ro',
        isa    => 'Int',
        writer => '_set_umask'
    );
    
    has 'uid' => (
        is     => 'ro',
        isa    => 'Int',
        writer => '_set_uid'
    );
    
    has 'log_dir' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_log_dir'
    );
    
    has 'log_file' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_log_file'
    );
    
    has '_log_file_in_use' => (
        is  => 'rw',
        isa => 'Bool'
    );
    
    has 'psgi_server' => (
        is      => 'ro',
        isa     => 'Object',
        lazy    => 1,
        builder => '_build_psgi_server'
    );
    
    has 'manager' => (
        is      => 'ro',
        isa     => 'VRPipe::Manager',
        lazy    => 1,
        builder => '_create_manager',        # we can't have a default because we must delay the create call
        handles => [qw(register_farm_server)]
    );
    
    has 'inmemory' => (
        is      => 'ro',
        isa     => 'VRPipe::Persistent::InMemory',
        lazy    => 1,
        builder => '_create_inmemory'
    );
    
    has '_warnings' => (
        is      => 'ro',
        isa     => 'ArrayRef',
        default => sub { [] },
        traits  => ['Array'],
        handles => {
            user_warning   => 'push',
            '_get_warning' => 'shift'
        }
    );
    
    has 'login_shell_script' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_login_shell_script'
    );
    
    method _build_schema {
        my $m = VRPipe::Manager->get;
        return $m->result_source->schema;
    }
    
    method _build_psgi_server {
        return Twiggy::Server::TLS->new(
            # host => `uname -n`, # parse the ip address, use as host? Currently defaults to 0.0.0.0
            port     => $self->port,
            tls_key  => $self->tls_key,
            tls_cert => $self->tls_cert
        );
    }
    
    method _create_manager {
        return VRPipe::Manager->create();
    }
    
    method _create_inmemory {
        return VRPipe::Persistent::InMemory->new();
    }
    
    sub BUILD {
        my $self = shift;
        
        # the idea is (per host) we'll have 1 daemon with 1 BackEnd instance per
        # dsn, and each dsn is supposed to have its own port associated with it.
        # we figure these out based on deployment
        my $deployment  = $self->deployment;
        my $vrp_config  = VRPipe::Config->new();
        my $method_name = $deployment . '_interface_port';
        my $port        = $vrp_config->$method_name();
        $port = "$port" if $port; # values retrieved from Config might be env vars, so we must force stringification
        unless ($port) {
            die "VRPipe SiteConfig had no port specified for $method_name\n";
        }
        $self->_set_port($port);
        
        # both production and testing will use the same tls key and cert
        my $tls_dir = $vrp_config->tls_dir;
        my $tls_key = file($tls_dir, $vrp_config->tls_key);
        $self->_set_tls_key("$tls_key");
        my $tls_cert = file($tls_dir, $vrp_config->tls_cert);
        $self->_set_tls_cert("$tls_cert");
        
        if (!-s $tls_key) {
            # create a private key file and certificate (we don't handle there
            # being no private key but a certificate; an existing cert will get
            # overwritten here!)
            system(qq[openssl req -new -x509 -sha256 -days 3650 -nodes -subj "/C=UK/O=VRPipe/OU=CN=VRPipe" -keyout $tls_key -out $tls_cert]);
            chmod(0600, $tls_key);
        }
        elsif (!-s $tls_cert) {
            # create a self-signed certificate using the existing key
            system(qq[openssl req -new -x509 -sha256 -days 3650 -key $tls_key -subj "/C=UK/O=VRPipe/OU=CN=VRPipe" -out $tls_cert]);
        }
        
        my $umask = $vrp_config->server_umask;
        $self->_set_umask("$umask");
        my $uid = $vrp_config->server_uid;
        $self->_set_uid("$uid");
        
        $method_name = $deployment . '_scheduler';
        my $scheduler = $vrp_config->$method_name();
        $self->_set_scheduler("$scheduler");
        
        my $login_shell_script = $vrp_config->login_shell_script();
        $self->_set_login_shell_script("$login_shell_script");
        
        VRPipe::Persistent::SchemaBase->database_deployment($deployment);
        $self->_set_dsn(VRPipe::Persistent::SchemaBase->get_dsn);
        
        $method_name = $deployment . '_logging_directory';
        my $log_dir      = $vrp_config->$method_name();
        my $log_basename = $self->dsn;
        $log_basename =~ s/\W/_/g;
        my $log_file = file($log_dir, 'vrpipe-server.' . $log_basename . '.log');
        $self->_set_log_file($log_file->stringify);
        $self->_set_log_dir("$log_dir");
        
        $method_name = $deployment . '_neo4j_server_url';
        my $neo4j_url = $vrp_config->$method_name();
        if ($neo4j_url) {
            $self->_set_neo4j_server_url("$neo4j_url");
        }
        
        require VRPipe::Persistent::Schema;
    }
    
    # mostly based on http://www.perlmonks.org/?node_id=374409, and used instead
    # of Proc::Daemon because we have slightly unusual requirements and have to
    # do things in our own special snow-flake way
    method daemonize {
        # -parent-
        
        # we fork and later exit the parent. That makes the launching process
        # think we're done and also makes sure we're not a process group leader
        my $pid = fork;
        if (defined $pid && $pid == 0) {
            # -first child-
            
            # we change directory to the root of all. That is a courtesy to the
            # system, which without it would be prevented from unmounting the
            # filesystem we started in. However, when testing we have to allow
            # for the possiblity of relative paths in fofn datasources, and we
            # have to know where t/ and modules/ are, so we assume that it is ok
            # to not chdir since the test will be over soon and the server auto-
            # stopped
            my $deployment = $self->deployment;
            unless ($deployment eq 'testing') {
                chdir '/' or die $!;
            }
            
            # we set the user-configured permissions mask for file creation
            umask $self->umask;
            
            # we call setsid() which does three things. We become leader of a
            # new session, group leader of a new process group, and become
            # detached from any terminal. We do the fork dance again to shuck
            # session leadership, which guarantees we never get a controlling
            # terminal.
            die "Cannot detach from controlling terminal" if setsid() < 0;
            $pid = fork;
            exit 0 if $pid;            # first child does not have to have to wait for the second child
            exit 1 if not defined $pid;
            
            # -second child-
            
            # we close unwanted filehandles we have inherited from the parent
            # process. That is another courtesy to the system, as well as being
            # sparing of our own process resources
            close $_ for (*STDIN, *STDOUT, *STDERR);
            
            # in production (when we might be running as root) set the real user
            # identifier and the effective user identifier for the daemon
            # process before opening files
            setuid($self->uid) if $deployment eq 'production';
            
            # reopen STDERR for logging purposes
            my $im = $self->inmemory();
            $im->log_stderr;
            
            return 1;
        }
        
        # -parent-
        
        # a child that terminates, but has not been waited for becomes a zombie.
        # So we wait for the first child to exit
        waitpid($pid, 0);
        exit 0;
    }
    
    method install_pipelines_and_steps {
        foreach my $module (findallmod(VRPipe::Pipelines)) {
            eval "require $module;";
            unless ($@) {
                my ($name) = $module =~ /VRPipe::Pipelines::(\w+)/;
                VRPipe::Pipeline->create(name => $name);
            }
        }
        
        foreach my $module (findallmod(VRPipe::Steps)) {
            eval "require $module;";
            unless ($@) {
                my ($name) = $module =~ /VRPipe::Steps::(\w+)/;
                VRPipe::Step->create(name => $name);
            }
        }
    }
    
    sub register_psgi_pages {
        my ($self, %pages) = @_;
        
        my $app = sub {
            my $env = shift; # PSGI env
            
            my $req      = Plack::Request->new($env);
            my $path     = $req->path_info;
            my @sub_args = ($env);
            my $sub      = $pages{$path};
            
            unless ($sub) {
                # see if the path matches any of the pages as a regex
                while (my ($regex, $this_sub) = each %pages) {
                    next unless $regex;
                    next if $regex eq '/';
                    if ($path =~ qr{$regex}) {
                        $sub = $this_sub;
                        push(@sub_args, "$1") if $1;
                        last;
                    }
                }
                keys %pages;
            }
            
            $sub ||= sub { $self->psgi_text_response(404, 'plain', '404: requested page (' . $req->request_uri . ') is not valid', $env); };
            
            my $return;
            eval { $return = &{$sub}(@sub_args); };
            if ($@) {
                my $err = $@;
                chomp($err);
                my $message = "fatal event captured responding to " . $req->request_uri . " for " . $req->address . ": " . $err;
                $self->log($message);
                my $res = $req->new_response(500);
                $res->content_type('text/plain');
                $res->body('500: ' . $message);
                return $res->finalize;
            }
            
            my $res;
            my $ref = ref($return);
            if ($ref) {
                if ($ref eq 'CODE') {
                    return $return;
                }
                $res = $return;
            }
            else {
                $res = $req->new_response(500);
                $res->content_type('text/plain');
                my $message = "page $path unexpectedly returned a non-reference";
                $self->log($message);
                $res->body('500: ' . $message);
            }
            
            return $res->finalize;
        };
        
        $self->psgi_server->register_service($app);
    }
    
    method psgi_text_response (Int $code, Str $format, Str $content, HashRef $env) {
        my $req = Plack::Request->new($env);
        my $res = $req->new_response($code);
        $res->content_type('text/' . $format);
        $res->body($content);
        return $res;
    }
    
    method psgi_file_response (Str $path, HashRef $env, ArrayRef :$regexes?, Int :$max_age?, Bool :$check_permissions = 0) {
        my $req = Plack::Request->new($env);
        
        # can I read it?
        my $ok = open(my $fh, $path);
        unless ($ok) {
            my $res = $req->new_response(404);
            $res->content_type('text/plain');
            if (!-e $path) {
                $res->body('404: ' . "$path not found");
            }
            else {
                $res->body('404: ' . "$path could not be opened: $!");
            }
            return $res;
        }
        
        if ($check_permissions) {
            # is the user allowed to read it?
            my $session = $req->cookies->{vrpipe_session};
            my $user;
            if ($session) {
                # check redis to see what session data we have for this user
                my $im           = $self->inmemory();
                my $session_hash = $im->get_session($session);
                
                # if we have data for the supplied session, and if there's a user,
                # they've previously authenticated successfully
                if ($session_hash && defined $session_hash->{user}) {
                    $user = $session_hash->{user};
                }
            }
            unless ($user) {
                my $res = $req->new_response(403);
                $res->content_type('text/plain');
                $res->body('403: For access to files you must first log in.');
                return $res;
            }
            
            # I can't just su as the user to see if they can read it, so have to
            # manually check the user and group and permissions on the file vs
            # the user. I deliberately do not check the parent directories since
            # I want to give access even when the user couldn't normally read
            # the file due to non-permissive directories
            my @stat     = stat($path);
            my $mode     = $stat[2];
            my $readable = 0;
            if ($mode & S_IROTH) {
                # it's world readable
                $readable = 1;
            }
            else {
                # see if the user owns the file
                my ($owner, undef, undef, $owner_gid) = getpwuid($stat[4]);
                if ($owner eq $user) {
                    # even if it's not strictly set to be user readable
                    # according to mode, we still allow access
                    $readable = 1;
                }
                elsif ($mode & S_IRGRP) {
                    # it's group readable; see if the user belongs in the file's
                    # group
                    my ($group, undef, $gid, $members) = getgrgid($stat[5]);
                    if ($owner_gid == $gid || $members =~ /\b$owner\b/) {
                        $readable = 1;
                    }
                }
            }
            
            unless ($readable) {
                my $res = $req->new_response(403);
                $res->content_type('text/plain');
                $res->body(q[403: You don't have permission to read this file.]);
                return $res;
            }
        }
        
        my $bin_mode = -B $path;
        my ($ext) = $path =~ /\.([^\.]+)$/;
        my $content_type;
        if ($ext && exists $extension_to_content_type{$ext}) {
            $content_type = $extension_to_content_type{$ext};
        }
        $content_type ||= $bin_mode ? $extension_to_content_type{bin} : $extension_to_content_type{txt};
        
        my $content = '';
        {
            local $/ = undef;
            
            if ($path =~ /vcf\.gz$/ && -s $path < 5000000) {
                close($fh);
                open($fh, "zcat $path |");
                $content_type = $extension_to_content_type{txt};
            }
            else {
                binmode $fh if $bin_mode;
            }
            $content = <$fh>;
            close($fh);
            
            foreach (@{ $regexes || [] }) {
                my ($search, $replace) = @{$_};
                $content =~ s/$search/$replace/g;
            }
        }
        
        my $res = $req->new_response(200);
        $res->content_type($content_type);
        $res->body($content);
        if ($max_age) {
            $res->header('cache-control' => "public, max-age=$max_age");
        }
        return $res;
    }
    
    sub psgi_nonblocking_json_response {
        my ($self, $graph, $sub, $env, @others) = @_;
        my $req = Plack::Request->new($env);
        
        return sub {
            my $responder = shift;
            
            fork_call {
                my $args;
                eval { $args = $graph->json_decode($req->content || '{}'); };
                if ($@) {
                    return $graph->json_encode({ errors => ['Unable to decode posted content'] });
                }
                if (ref($args) ne 'HASH') {
                    return $graph->json_encode({ errors => ['Posted content was not a hash'] });
                }
                
                # I couldn't get Plack::Middleware::Session to work, and want
                # the sessions in redis anyway, so I implement my own session
                # management by passing my own vrpipe_session cookie value if
                # it exists, then when we return from this fork I check to see
                # if the returned json had vrpipe_session and set it as a
                # cookie if so
                $args->{vrpipe_session} = $req->cookies->{vrpipe_session};
                
                my $user;
                my $session = $args->{vrpipe_session};
                if ($session) {
                    my $im = $self->inmemory();
                    
                    # check redis to see what session data we have for this user
                    my $session_hash = $im->get_session($session);
                    
                    # if we have data for the supplied session, and if there's a user,
                    # they've previously authenticated successfully
                    if ($session_hash && defined $session_hash->{user}) {
                        $user = $session_hash->{user};
                    }
                    else {
                        # their session expired in some way
                        delete $args->{vrpipe_session};
                    }
                }
                $args->{authenticated_user} = $user;
                
                my $data;
                eval { $data = &{$sub}($args, @others); };
                if ($@) {
                    $data = { errors => [$@] };
                }
                
                return $graph->json_encode($data);
            }
            sub {
                my ($json) = @_;
                my $res = $req->new_response(200);
                $res->content_type('application/json');
                
                my $decoded = $graph->json_decode($json);
                my $session = delete $decoded->{vrpipe_session} if ref($decoded) eq 'HASH';
                if (defined $session) {
                    $res->cookies->{vrpipe_session} = { value => $session, path => "/" };
                    $json = $graph->json_encode($decoded);
                }
                
                $res->body($json);
                
                $responder->($res->finalize);
            };
        };
    }
    
    method req_to_opts (Object $req, ArrayRef $ctp_array?) {
        my %opts = %{ $req->parameters->mixed || {} };
        
        my $convert_to_persistent = delete $opts{convert_to_persistent} || $ctp_array;
        if ($convert_to_persistent) {
            unless (ref($convert_to_persistent)) {
                $convert_to_persistent = [$convert_to_persistent];
            }
            foreach my $ctp (@$convert_to_persistent) {
                my ($opt, $class) = split(/\!/, $ctp);
                next unless defined $opts{$opt};
                my $val        = $opts{$opt};
                my @desired    = ref($val) eq 'ARRAY' ? @$val : ($val);
                my $full_class = "VRPipe::$class";
                
                my @found;
                foreach my $desired (@desired) {
                    my $found;
                    if ($desired =~ /^\d+$/) {
                        ($found) = $full_class->search({ id => $desired });
                        unless ($found) {
                            die "$desired is not a valid $class id\n";
                        }
                    }
                    else {
                        ($found) = $full_class->search({ name => $desired });
                        unless ($found) {
                            die "$desired is not a valid $class name\n";
                        }
                    }
                    push(@found, $found);
                }
                
                if (ref($val) eq 'ARRAY') {
                    $opts{$opt} = \@found;
                }
                else {
                    $opts{$opt} = $found[0];
                    if (@found > 1) {
                        $self->user_warning("--$opt @desired resulted in more than one $class object; only using the first");
                    }
                }
            }
        }
        
        foreach my $key (keys %opts) {
            if ($key =~ /^(.+)\:hashderef\:(.+)$/) {
                my $val = delete $opts{$key};
                $opts{$1}->{$2} = $val;
            }
        }
        
        return \%opts;
    }
    
    method get_pipelinesetups (HashRef $opts, Bool $inactive?, Bool $allow_no_setups?) {
        my $multi_setups = $opts->{'_multiple_setups'};
        my @requested_setups = defined $opts->{setup} ? ($multi_setups ? @{ ref($opts->{setup}) eq 'ARRAY' ? $opts->{setup} : [$opts->{setup}] } : ($opts->{setup})) : ();
        
        my @setups;
        if (@requested_setups) {
            @setups = @requested_setups;
        }
        elsif ($multi_setups) {
            my $user = $opts->{user} || 'all';
            unless (defined $inactive) {
                $inactive = $opts->{deactivated};
            }
            @setups = VRPipe::PipelineSetup->search({ $user eq 'all' ? () : (user => $user), $inactive ? () : (active => 1) }, { prefetch => ['datasource', 'pipeline'] });
        }
        
        if ($multi_setups && !@setups && !$allow_no_setups) {
            die "No PipelineSetups match your settings (did you remember to specify --user?)\n";
        }
        
        return $multi_setups ? @setups : $setups[0];
    }
    
    method sub_modules (Str $base) {
        $base = "VRPipe::$base";
        my @modules = findsubmod $base;
        unless (@modules) {
            @modules = findsubmod "${base}s";
        }
        my @names;
        foreach my $module (@modules) {
            my ($name) = $module =~ /::([^:]+)$/;
            push(@names, $name);
        }
        return @names;
    }
    
    method make_all_objects (Str $class) {
        my @modules = $self->sub_modules($class);
        $class = "VRPipe::$class";
        foreach my $name (@modules) {
            $class->get(name => $name);
        }
    }
    
    method ask_for_object (Str :$question!, Str :$class!, Str :$column!) {
        $self->make_all_objects($class);
        my @things     = "VRPipe::$class"->search({});
        my %things     = map { $_->$column => $_ } @things;
        my @thing_keys = sort keys %things;
        $self->output("\n");
        my %num_to_key;
        foreach my $i (0 .. $#thing_keys) {
            my $num = $i + 1;
            my $key = $thing_keys[$i];
            $num_to_key{$num} = $key;
            my $output = "$num. $key";
            my $obj    = $things{$key};
            if ($obj->can('description')) {
                $output .= ' [' . $obj->description . ']';
            }
            $self->output($output);
        }
        my $chosen_num = $self->pick_number(question => $question, max => scalar(@thing_keys));
        return $things{ $num_to_key{$chosen_num} };
    }
    
    method already_exists (Str $class!, Str $key!, Str $value!) {
        my $found = "VRPipe::$class"->search({ $key => $value });
        if ($found) {
            return "a $class already exists with $key '$value'";
        }
        return;
    }
    
    # we have our own ssh wrapper (instead of using Net::SSH) because we want
    # to supply ssh options and handle stderr/out and error handling ourselves.
    # NB: $cmd will be surrounded by double quotes, so cannot contain double
    # quotes or unescaped $variables
    method ssh (Str $host, Str $cmd, Str :$working_dir?) {
        my $background = !defined wantarray();
        my $ssh_opts   = '-T -n -o BatchMode=yes -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -o LogLevel=quiet -o ConnectionAttempts=1 -o ConnectTimeout=5';
        my $ssh_cmd    = qq[ssh $ssh_opts $host ];
        $ssh_cmd .= q["];
        
        # we must be sure to get the user's expected environment variables, so
        # we'll source their shell login script (which doesn't happen in our
        # non-interactive tty-less ssh invocation)
        my $login_shell_script = $self->login_shell_script;
        if ($login_shell_script && -s $login_shell_script) {
            $ssh_cmd .= qq[source $login_shell_script; ];
        }
        
        # we might need to cd do a directory before executing the command
        if ($working_dir) {
            $ssh_cmd .= qq[cd $working_dir; ];
        }
        
        $ssh_cmd .= q[nohup ] if $background;
        $ssh_cmd .= $cmd;
        $ssh_cmd .= q[ &>/dev/null] if ($background && $cmd !~ /&>/);
        $ssh_cmd .= q[ &]           if $background;
        $ssh_cmd .= q["];
        
        unless ($background) {
            #*** we could do it like Net::SSH does it with an open3 call to
            # capture the stderr, but...
            # my $tied = tied *STDERR ? 1 : 0;
            # untie *STDERR if $tied;
            # my $return;
            # ...
            # $self->log_stderr() if $tied;
            
            # ignore errors (stderr will just end up in logs) and assume the
            # caller will check the return value is sensible
            return `$ssh_cmd`;
        }
        else {
            my $system = system($ssh_cmd);
            unless ($system == 0) {
                $self->log("ssh failed: $?");
            }
            return;
        }
    }
}

1;
