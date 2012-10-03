
=head1 NAME

VRPipe::Interface::CmdLine - shared methods of all cmd-line frontend scripts

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

The various B<vrpipe-*> scripts use this module which provides methods that
they all need to function.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Interface::CmdLine {
    use Getopt::Long qw(GetOptions GetOptionsFromString);
    use Perl6::Form;
    use LWP::UserAgent;
    use Sys::Hostname::Long;
    use VRPipe::Config;
    use VRPipe::Persistent::SchemaBase;
    use Config;
    
    has 'description' => (
        is       => 'rw',
        isa      => 'Str',
        required => 1
    );
    
    has 'extra_args' => (
        is  => 'rw',
        isa => 'Str'
    );
    
    has 'opt_spec' => (
        is      => 'rw',
        isa     => 'ArrayRef[ArrayRef]',
        lazy    => 1,
        builder => '_default_opt_spec'
    );
    
    has '_opts_hash' => (
        is      => 'ro',
        isa     => 'HashRef',
        default => sub { {} },
        writer  => '_set_opts',
        traits  => ['Hash'],
        handles => {
            opts           => 'get',
            option_was_set => 'defined',
            '_set_opt'     => 'set'
        }
    );
    
    has 'usage' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_usage'
    );
    
    has '_multiple_setups' => (
        is  => 'rw',
        isa => 'Bool'
    );
    
    has 'no_user_option' => (
        is      => 'ro',
        isa     => 'Bool',
        default => 0
    );
    
    has '_ua' => (
        is      => 'ro',
        isa     => 'LWP::UserAgent',
        lazy    => 1,
        builder => '_build_ua'
    );
    
    has 'port' => (
        is     => 'ro',
        isa    => PositiveInt,
        writer => '_set_port'
    );
    
    has 'dsn' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_dsn'
    );
    
    has 'server_ok' => (
        is      => 'ro',
        isa     => 'Bool',
        lazy    => 1,
        builder => 'check_server'
    );
    
    has '_ua_port_baseurl' => (
        is      => 'ro',
        isa     => 'ArrayRef',
        lazy    => 1,
        builder => '_build_ua_port_baseurl'
    );
    
    method _default_opt_spec {
        return [['deployment=s', 'Use the production or testing database', { default => 'production' }], ['env|e=s', 'Use options stored in an environment variable'], ['help|h', 'Print this usage message and exit']];
    }
    
    method _build_ua {
        my $ua = LWP::UserAgent->new(timeout => 20, agent => 'VRPipe-Client');
        return $ua;
    }
    
    sub BUILD {
        my $self = shift;
        
        # we initially used Getopt::Long::Descriptive, but because it does not
        # show what kind of arg each option takes, we implement here a similar
        # interface, but do not use the actual GLD code.
        
        my $opt_spec = $self->opt_spec;
        unless (@$opt_spec <= 4 && @$opt_spec >= 3 && $opt_spec->[-1]->[0] eq 'help|h') {
            my $default = $self->_default_opt_spec;
            
            foreach my $opt_spec (@$opt_spec) {
                my ($def, $help, $extra) = @$opt_spec;
                if ($help) {
                    my ($name, $req_or_opt, $type) = split(/([=:])/, $def);
                    if ($name eq 'setup') {
                        if ($type eq 's@') {
                            unless ($self->no_user_option) {
                                $default->[4] = $default->[2];
                                $default->[2] = ['user|u=s', 'Only show entries for PipelineSetups created by this user; use "all" to show entries for all users', { default => getlogin || getpwuid($<) || 'vrpipe' }];
                                $default->[3] = ['deactivated', 'Also show deactivated PipelineSetups', { default => 0 }];
                            }
                            $self->_multiple_setups(1);
                        }
                    }
                }
            }
            push(@$opt_spec, [], ['General options:'], @$default);
        }
        
        my @opts;
        my %defaults;
        my %convert_to_persistent;
        my $script_name = file($0)->basename;
        my $usage       = '';
        my (@shorts, $has_long);
        foreach my $opt_spec (@{ $self->opt_spec }) {
            my ($def, $help, $extra) = @$opt_spec;
            if ($help) {
                push(@opts, $def);
                
                my ($name, $req_or_opt, $type) = split(/([=:])/, $def);
                
                my ($long, $short) = length($name) == 1 ? (undef, $name) : (split(/\|/, $name));
                if (!$short && length($long) == 1) {
                    $short = $long;
                    undef($long);
                }
                if ($short && length($short) != 1) {
                    undef($short);
                }
                my $option = '';
                if ($short) {
                    push(@shorts, $short);
                    $option .= ' -' . $short;
                }
                if ($long) {
                    $has_long = 1;
                    $option .= ' --' . $long;
                }
                
                if ($extra) {
                    my $default = $extra->{default};
                    if (defined $default) {
                        $defaults{ $long || $short } = $default;
                        $help = "[$default] " . $help;
                    }
                    
                    my $persistent_class = $extra->{persistent_object};
                    if ($persistent_class) {
                        $convert_to_persistent{ $long || $short } = $persistent_class;
                    }
                    
                    next if $extra->{hidden};
                }
                
                my $value = '';
                if ($type) {
                    if ($type =~ /^s/) {
                        $value = '<str>';
                    }
                    elsif ($type =~ /^i/) {
                        $value = '<int>';
                    }
                    elsif ($type =~ /^f/) {
                        $value = '<num>';
                    }
                    $value .= $req_or_opt eq ':' ? '?' : '';
                }
                
                $usage .= form "  {[[[[[[[[[[[[[[[} {IIII} {[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[}", $option, $value, $help;
            }
            elsif ($def) {
                $usage .= form "{[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[}", $def;
            }
            else {
                $usage .= "\n";
            }
        }
        my $o = '';
        if (@shorts) {
            $o .= ' [-' . join('', sort @shorts) . ']';
        }
        if ($has_long) {
            $o .= ' [long options...]';
        }
        my $extra_args = $self->extra_args;
        $o .= ' ' . $extra_args if $extra_args;
        $usage = $self->description . "\n$script_name$o\n" . $usage;
        $self->_set_usage($usage);
        
        my %opts;
        Getopt::Long::Configure("bundling");
        $self->help unless GetOptions(\%opts, @opts);
        $self->_set_opts(\%opts);
        
        my $from_env = $self->opts('env');
        if ($from_env) {
            $self->help unless GetOptionsFromString($ENV{$from_env}, \%opts, @opts);
            $self->_set_opts(\%opts);
        }
        
        while (my ($opt, $val) = each %defaults) {
            next if $self->option_was_set($opt);
            $self->_set_opt($opt => $val);
        }
        
        my @ctps;
        while (my ($opt, $class) = each %convert_to_persistent) {
            push(@ctps, "$opt!$class");
        }
        if (@ctps) {
            $self->_set_opt(convert_to_persistent => \@ctps);
        }
        
        my $help       = $self->opts('help');
        my $deployment = $self->opts('deployment');
        unless ($deployment eq 'production' || $deployment eq 'testing') {
            $self->error("--deployment must be <production|testing>");
            $help = 1;
        }
        $self->help if $help;
        
        my $vrp_config  = VRPipe::Config->new();
        my $method_name = $deployment . '_interface_port';
        my $port        = $vrp_config->$method_name();
        unless ($port) {
            $self->throw("VRPipe SiteConfig had no port specified for $method_name");
        }
        $self->_set_port("$port");                                        # values retrieved from Config might be env vars, so we must force stringification
        VRPipe::Persistent::SchemaBase->database_deployment($deployment); # this does not access the db, just lets get_dsn method work
        $self->_set_dsn(VRPipe::Persistent::SchemaBase->get_dsn);
    }
    
    method help {
        $self->die_with_error($self->usage);
    }
    
    method output (@messages) {
        chomp($messages[-1]);
        print @messages, "\n";
    }
    
    method error (@messages) {
        chomp($messages[-1]);
        warn @messages, "\n";
    }
    
    method die_with_error (@messages) {
        $self->error(@messages);
        die "\n";
    }
    
    method _build_ua_port_baseurl {
        my $ua       = $self->_ua;
        my $port     = $self->port;
        my $base_url = 'http://' . hostname_long . ':' . $port;
        return [$ua, $port, $base_url];
    }
    
    method check_server (Bool $no_auto_start_or_die = 0) {
        my ($ua, $port, $base_url) = @{ $self->_ua_port_baseurl };
        
        # try and get a response from the port
        my @post_args = ($base_url . '/dsn');
        my $response  = $ua->post(@post_args);
        my $server_dsn;
        if ($response->is_success) {
            $server_dsn = $response->decoded_content;
        }
        elsif ($no_auto_start_or_die) {
            return 0;
        }
        else {
            if ($response->code == 500) {
                $self->error("Can't connect to VRPiper server at $base_url, will attempt to auto-start it...");
                my $script = $self->vrpipe_script_path('vrpipe-server', $self->opts('deployment'));
                my $cmd = $script . ' start'; # (this confuses perltidy if put directly in the system() call)
                system($cmd);
                
                # the vrpipe-server call returns ~instantly, but may take some
                # time before the server is actually ready to connect to; keep
                # trying for the next 20 seconds
                my $seconds = 20;
                while ($seconds--) {
                    $response = $ua->post(@post_args);
                    if ($response->is_success) {
                        $server_dsn = $response->decoded_content;
                        last;
                    }
                    sleep(1);
                }
            }
            $self->die_with_error($response->status_line) unless $server_dsn;
        }
        
        # make sure the response is valid
        my $expected_dsn = $self->dsn;
        unless ($server_dsn eq $expected_dsn) {
            if ($no_auto_start_or_die) {
                $self->error("expected dsn $expected_dsn but got $server_dsn");
                return -1;
            }
            else {
                $self->die_with_error("There is a server bound to port $port, but it is either not connected to the correct database ($expected_dsn), or is not a VRPipe server at all.\nIts reported dsn was:\n$server_dsn\n");
            }
        }
        
        return 1;
    }
    
    method server_get (Str $page, HashRef $opts?) {
        $opts ||= {};
        $opts->{display_format} = 'plain';
        
        # first we need to make sure that the server bound to the port we
        # configured for our database is actually connected to our database
        # (and we'll auto-start the server if it isn't running at all)
        $self->die_with_error("Can't connect to server.") unless $self->server_ok;
        
        # now we'll handle the desired request
        my ($ua, undef, $base_url) = @{ $self->_ua_port_baseurl };
        my @post_args = ($base_url . $page, [$self->_deref_hash($self->_opts_hash), '_multiple_setups' => $self->_multiple_setups, $self->_deref_hash($opts)]);
        
        my $response = $ua->post(@post_args);
        if ($response->is_success) {
            return $response->decoded_content;
        }
        else {
            $self->die_with_error("Failed to get $post_args[0] (@{$post_args[1]}):" . $response->status_line);
        }
    }
    
    method _deref_hash (HashRef $hash) {
        my @list;
        while (my ($key, $val) = each %$hash) {
            my $ref = ref($val);
            if ($ref) {
                if ($ref eq 'ARRAY') {
                    foreach my $subval (@$val) {
                        push(@list, $key, $subval);
                    }
                }
                elsif ($ref eq 'HASH') {
                    while (my ($sub_key, $subval) = each %{$val}) {
                        push(@list, "$key:hashderef:$sub_key", $subval);
                    }
                }
                else {
                    $self->throw("Can't cope with values of type $ref as post args");
                }
            }
            else {
                push(@list, $key, $val);
            }
        }
        return @list;
    }
    
    method ask_question (Str :$question!, ArrayRef :$possibles?, Str :$allow_multiple?, Str :$default?, Bool :$required?, CodeRef :$not_allowed?, ArrayRef :$na_args?) {
        undef $possibles unless $possibles && @$possibles;
        if (defined $default && length($default) == 0) {
            undef $default;
        }
        
        print STDERR "$question";
        my %allowed;
        if ($possibles) {
            print STDERR " <", join('|', @$possibles), ">";
            %allowed = map { $_ => 1 } @$possibles;
        }
        if (defined $default) {
            print STDERR " [$default]";
        }
        print STDERR ": ";
        
        my $answer;
        do {
            $answer = <STDIN>;
            chomp($answer);
            
            if ($possibles) {
                my $valid = 1;
                my @answers = $allow_multiple ? split(/$allow_multiple/, $answer) : ($answer);
                foreach my $sub_answer (@answers) {
                    unless (exists $allowed{$sub_answer}) {
                        $valid = 0;
                    }
                }
                
                unless ($valid) {
                    if (defined $default && !$answer) {
                        $answer = $default;
                    }
                    else {
                        undef $answer;
                        my $one_of = $allow_multiple ? 'one or more of' : 'one of';
                        print STDERR "Your answer can only consist of $one_of <", join('|', @$possibles), ">: ";
                    }
                }
            }
            elsif (!$answer && "$answer" ne "0") {
                if (defined $default) {
                    $answer = $default;
                }
                elsif ($required) {
                    undef $answer;
                    print STDERR "An answer is required. Try again: ";
                }
            }
            
            if ($answer && $not_allowed) {
                push(@$na_args, $answer);
                my $reason = &$not_allowed(@$na_args);
                if ($reason) {
                    undef $answer;
                    pop(@$na_args);
                    print STDERR "Your answer isn't allowed because $reason. Try again: ";
                }
            }
        } while (!defined $answer);
        
        return $answer;
    }
    
    method pick_number (Str :$question!, PositiveInt :$max!, PositiveInt :$default?) {
        return $self->ask_question(question => $question, possibles => [1 .. $max], required => 1, $default ? (default => $default) : ());
    }
    
    # we might not have vrpipe-* in our PATH, and the modules it needs might not
    # be in our PERL5LIB, so allow us to still work if we're running from the
    # git repo root dir (eg. during testing prior to an install). In fact, for
    # testing purposes, we prefer this to some version of the files installed
    # elsewhere.
    method vrpipe_script_path (Str $script, Str $deployment) {
        my $path;
        if ($deployment eq 'testing') {
            my $local_script = File::Spec->catfile('scripts', $script);
            if (-x $local_script && -d 'modules' && -d 't') {
                my $thisperl = $Config{perlpath};
                if ($^O ne 'VMS') {
                    $thisperl .= $Config{_exe} unless $thisperl =~ m/$Config{_exe}$/i;
                }
                $path = "$thisperl -Imodules -It $local_script";
            }
        }
        else {
            $path = $script;
        }
        
        if ($deployment eq 'testing') {
            $path .= " --deployment testing";
        }
        
        return $path;
    }
}

1;
