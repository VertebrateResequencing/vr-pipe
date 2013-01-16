
=head1 NAME

VRPipe::Interface::BackEnd - methods for the server to handle user interaction

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

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

class VRPipe::Interface::BackEnd {
    use Perl6::Form;
    use Module::Find;
    use VRPipe::Persistent::SchemaBase;
    use VRPipe::Config;
    use AnyEvent::ForkManager;
    use Twiggy::Server;
    use Continuity;
    use Continuity::Adapt::PSGI;
    use Plack::Request;
    use XML::LibXSLT;
    use XML::LibXML;
    use Sys::Hostname::Long;
    use POSIX qw(setsid setuid);
    use Cwd qw(chdir getcwd);
    use Time::Format;
    use Module::Find;
    use Email::Sender::Simple;
    use Email::Simple::Creator;
    use Fcntl qw/:flock SEEK_END/;
    
    my $xsl_html = <<'XSL';
<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="html"/>

<xsl:template match="/">
    <html>
        <head>
            <title><xsl:value-of select="/interface/title"/></title>
            
            <style type="text/css">
                dl {
                    font-size: 9pt
                }
                dl dt {
                    background:#99CCFF;
                    color:#fff;
                    float:left;
                    font-weight:bold;
                    margin-right:10px;
                    padding:5px;
                    width:65px;
                }
                dl dd {
                    margin:2px 0;
                    padding:5px 0;
                }
            </style>
        </head>
        
        <body>
            <xsl:apply-templates/>  
        </body>
    </html>
</xsl:template>

<xsl:template match="/interface/error">
    <p style="color:red"><xsl:value-of select="."/></p>
</xsl:template>

<xsl:template match="/interface/warning"/>

<xsl:template match="/interface/title">
    <h1><xsl:value-of select="."/></h1>
    
    <xsl:for-each select="/interface/warning">
        <p><span style="background-color:yellow"><xsl:value-of select="."/></span></p>
    </xsl:for-each>
</xsl:template>

<xsl:template match="/interface/response_line">
    <p><xsl:value-of select="."/></p>
</xsl:template>

<xsl:template match="objects">
    <table border="1">
        <thead>
            <tr bgcolor="#9acd32">
                <xsl:for-each select="./object[1]/attribute">
                    <th><xsl:value-of select="@name"/></th>
                </xsl:for-each>
            </tr>
        </thead>
        <tbody>
            <xsl:apply-templates select="object"/>
        </tbody>
    </table>
</xsl:template>

<xsl:template match="object">
    <tr>
        <xsl:apply-templates select="attribute"/>
    </tr>
</xsl:template>

<xsl:template match="attribute[object]">
    <td>
        <dl class="object_attribs">
            <xsl:for-each select="./object/attribute">
                <dt><xsl:value-of select="@name"/></dt>
                <dd>
                    <xsl:choose>
                        <xsl:when test="./object">
                            <xsl:value-of select="./object/attribute[@name='id']"/>
                        </xsl:when>
                        <xsl:when test="./hash">
                            <xsl:for-each select="./hash/pair">
                                <b><xsl:value-of select="@key"/>: </b><xsl:value-of select="."/><br/>
                            </xsl:for-each>
                        </xsl:when>
                        <xsl:otherwise>
                            <xsl:value-of select="."/>
                        </xsl:otherwise>
                    </xsl:choose>
                </dd>
            </xsl:for-each>
        </dl>
    </td>
</xsl:template>

<xsl:template match="attribute[hash]">
    <td>
        <xsl:apply-templates select="hash"/>
    </td>
</xsl:template>

<xsl:template match="attribute">
    <td>
        <xsl:value-of select="."/>
    </td>
</xsl:template>

<xsl:template match="hash">
    <dl class="hash">
        <xsl:apply-templates select="pair"/>
    </dl>
</xsl:template>

<xsl:template match="pair">
    <dt><xsl:value-of select="@key"/></dt>
    <dd><xsl:value-of select="."/></dd>
</xsl:template>

</xsl:stylesheet>
XSL
    
    my $xsl_plain = <<'XSL';
<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text"/>

<xsl:template match="/">
    <xsl:apply-templates/>
</xsl:template>

<xsl:template match="/interface/error">
    <xsl:value-of select="."/><xsl:text>
</xsl:text>
</xsl:template>

<xsl:template match="/interface/title"/>

<xsl:template match="/interface/warning">
    <xsl:value-of select="."/><xsl:text>
</xsl:text>
</xsl:template>

<xsl:template match="/interface/response_line">
    <xsl:value-of select="."/><xsl:text>
</xsl:text>
</xsl:template>

<xsl:template match="objects">
    <xsl:apply-templates select="object"/>
</xsl:template>

<xsl:template match="object[@class='PipelineSetup']">
    <xsl:text> --- Pipeline Setup '</xsl:text>
    <xsl:value-of select="./attribute[@name='name']"/>
    <xsl:text>' (id </xsl:text>
    <xsl:value-of select="./attribute[@name='id']"/>
    <xsl:text> for user </xsl:text>
    <xsl:value-of select="./attribute[@name='user']"/>
    <xsl:text>)</xsl:text>
    <xsl:if test="./attribute[@name='active'] = 0">
        <xsl:text> currently DEACTIVATED</xsl:text>
    </xsl:if>
    <xsl:text> ---
</xsl:text>
    
    <xsl:if test="(@display_mode!='list') and (@display_mode!='defunct')">
        <xsl:if test="@display_mode='full'">
            <xsl:text>Pipeline: </xsl:text>
            <xsl:value-of select="./attribute[@name='pipeline']/object/attribute[@name='name']"/>
            <xsl:text> | </xsl:text>
            <xsl:value-of select="./attribute[@name='pipeline']/object/attribute[@name='num_steps']"/>
            <xsl:text> steps | </xsl:text>
            <xsl:value-of select="./attribute[@name='pipeline']/object/attribute[@name='description']"/>
            <xsl:text>
</xsl:text>
            <xsl:text>PipelineSetup options:
</xsl:text>
            <xsl:choose>
                <xsl:when test="./attribute[@name='options']/hash">
                    <xsl:apply-templates select="./attribute[@name='options']/hash/pair"/>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:value-of select="./attribute[@name='options']"/>
                    <xsl:text>
</xsl:text>
                </xsl:otherwise>
            </xsl:choose>
            <xsl:text>PipelineSetup output root: </xsl:text>
            <xsl:value-of select="./attribute[@name='output_root']"/>
            <xsl:text>
</xsl:text>
            <xsl:text>DataSource: </xsl:text>
            <xsl:value-of select="./attribute[@name='datasource']/object/attribute[@name='id']"/>
            <xsl:text> | </xsl:text>
            <xsl:value-of select="./attribute[@name='datasource']/object/attribute[@name='type']"/>
            <xsl:text> | </xsl:text>
            <xsl:value-of select="./attribute[@name='datasource']/object/attribute[@name='method']"/>
            <xsl:text> | </xsl:text>
            <xsl:value-of select="./attribute[@name='datasource']/object/attribute[@name='source']"/>
            <xsl:text>
</xsl:text>
            <xsl:choose>
                <xsl:when test="./attribute[@name='datasource']/object/attribute[@name='options']/hash">
                    <xsl:apply-templates select="./attribute[@name='datasource']/object/attribute[@name='options']/hash/pair"/>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:value-of select="./attribute[@name='datasource']/object/attribute[@name='options']"/>
                    <xsl:text>
</xsl:text>
                </xsl:otherwise>
            </xsl:choose>
            
            <xsl:text>
</xsl:text>
        </xsl:if>
        
        <xsl:text>There are a total of </xsl:text>
        <xsl:value-of select="./attribute[@name='elements_total']"/>
        <xsl:text> Data Elements in the datasource to work on, and </xsl:text>
        <xsl:value-of select="./attribute[@name='elements_incomplete']"/>
        <xsl:text> elements are incomplete
</xsl:text>
        <xsl:choose>
            <xsl:when test="./attribute[@name='steps_completed']">
                <xsl:text>Breakdown:
</xsl:text>
                <xsl:for-each select="./attribute[@name='steps_completed']/hash/pair">
                    <xsl:text>    </xsl:text>
                    <xsl:value-of select="@key"/>
                    <xsl:text> steps completed => </xsl:text>
                    <xsl:value-of select="."/>
                    <xsl:text>
</xsl:text>
                </xsl:for-each>
                <xsl:value-of select="./attribute[@name='completion']/@explanation"/>
                
                <xsl:if test="./attribute[@name='submission_state']">
                    <xsl:text>
                    
Current submission state:
</xsl:text>
                    <xsl:apply-templates select="./attribute[@name='submission_state']/hash/pair"/>
                </xsl:if>
            </xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="./attribute[@name='completion']/@explanation"/>
                <xsl:text>
</xsl:text>
            </xsl:otherwise>
        </xsl:choose>
        
        <xsl:if test="./attribute[@name='problems']">
            <xsl:value-of select="./attribute[@name='problems']"/>
            <xsl:text>
</xsl:text>
        </xsl:if>
        
        <xsl:text>------

</xsl:text>
    </xsl:if>
    <xsl:if test="@display_mode='defunct'">
        <xsl:if test="./attribute[@name='problems']">
            <xsl:value-of select="./attribute[@name='problems']"/>
            <xsl:text>
</xsl:text>
        </xsl:if>
    </xsl:if>
</xsl:template>

<xsl:template match="pair">
    <xsl:text>    </xsl:text>
    <xsl:value-of select="@key"/>
    <xsl:text> => </xsl:text>
    <xsl:value-of select="."/>
    <xsl:text>
</xsl:text>
</xsl:template>

</xsl:stylesheet>
XSL
    
    my $parser = XML::LibXML->new();
    my $xslt   = XML::LibXSLT->new();
    our %display_format_to_stylesheet = (
        html  => $xslt->parse_stylesheet($parser->load_xml(string => $xsl_html)),
        plain => $xslt->parse_stylesheet($parser->load_xml(string => $xsl_plain))
    );
    
    has 'deployment' => (
        is       => 'ro',
        isa      => 'Str',
        required => 1
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
    
    has 'scheduler' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_scheduler'
    );
    
    has 'admin_user' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_admin_user'
    );
    
    has 'email_domain' => (
        is     => 'ro',
        isa    => 'Str',
        writer => '_set_email_domain'
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
    
    has 'continuations' => (
        is      => 'ro',
        isa     => 'HashRef',
        default => sub { {} },
        traits  => ['Hash'],
        handles => {
            store_continuation       => 'set',
            continuation_for_page    => 'get',
            continuation_page_stored => 'exists'
        }
    );
    
    has 'manager' => (
        is      => 'ro',
        isa     => 'VRPipe::Manager',
        lazy    => 1,
        builder => '_create_manager',        # we can't have a default because we must delay the create call
        handles => [qw(register_farm_server)]
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
    
    method _build_schema {
        my $m = VRPipe::Manager->get;
        return $m->result_source->schema;
    }
    
    method _build_psgi_server {
        my $port = $self->port;
        return Twiggy::Server->new(port => $port); # host `uname -n`; parse the ip address, use as host? Currently defaults to 0.0.0.0
    }
    
    method _create_manager {
        return VRPipe::Manager->create();
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
        unless ($port) {
            die "VRPipe SiteConfig had no port specified for $method_name\n";
        }
        $self->_set_port("$port"); # values retrieved from Config might be env vars, so we must force stringification
        
        my $umask = $vrp_config->server_umask;
        $self->_set_umask("$umask");
        my $uid = $vrp_config->server_uid;
        $self->_set_uid("$uid");
        
        $method_name = $deployment . '_scheduler';
        my $scheduler = $vrp_config->$method_name();
        $self->_set_scheduler("$scheduler");
        
        my $admin_user = $vrp_config->admin_user();
        $self->_set_admin_user("$admin_user");
        
        my $email_domain = $vrp_config->email_domain();
        $self->_set_email_domain("$email_domain");
        
        VRPipe::Persistent::SchemaBase->database_deployment($deployment);
        $self->_set_dsn(VRPipe::Persistent::SchemaBase->get_dsn);
        
        $method_name = $deployment . '_logging_directory';
        my $log_dir      = $vrp_config->$method_name();
        my $log_basename = $self->dsn;
        $log_basename =~ s/\W/_/g;
        my $log_file = file($log_dir, 'vrpipe-server.' . $log_basename . '.log');
        $self->_set_log_file($log_file->stringify);
        
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
            $self->log_stderr;
            
            return 1;
        }
        
        # -parent-
        
        # a child that terminates, but has not been waited for becomes a zombie.
        # So we wait for the first child to exit
        waitpid($pid, 0);
        exit 0;
    }
    
    method log_stderr {
        my $ok = open(STDERR, '>>', $self->log_file);
        if ($ok) {
            $self->_log_file_in_use(1);
        }
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
            
            my $req  = Plack::Request->new($env);
            my $page = $req->path_info;
            my $sub  = $pages{$page} || sub { $self->psgi_text_response(404, 'plain', '404: requested page (' . $req->request_uri . ') is not valid', $env); };
            
            my $return;
            eval { $return = &{$sub}($env); };
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
                my $message = "page $page unexpectedly returned a non-reference";
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
    
    sub psgi_nonblocking_xml_response {
        my ($self, $sub, $env, @others) = @_;
        my $req = Plack::Request->new($env);
        
        return sub {
            my $responder = shift;
            
            # (we use AnyEvent::ForkManager instead of AnyEvent::Util::fork_call
            #  because the latter doesn't seem to actually fork or do anything
            #  async)
            my $pm = AnyEvent::ForkManager->new(max_workers => 1);
            
            $pm->start(
                cb => sub {
                    my $xml;
                    try {
                        $xml = &{$sub}($req, @others);
                    }
                    catch ($err) {
                        chomp($err);
                        $xml = $self->xml_tag('error', $err);
                        $self->log("fatal event captured responding to " . $req->request_uri . " for " . $req->address . ": " . $err);
                    }
                    
                    my $warnings = '';
                    while (my $warning = $self->_get_warning) {
                        $warnings .= $self->xml_tag('warning', $warning);
                    }
                    
                    my $format = $req->param('display_format');
                    $format ||= 'html';
                    my $content = $self->transform_xml($warnings . $xml, $format);
                    
                    my $res = $req->new_response(200);
                    $res->content_type('text/' . $format);
                    $res->body($content);
                    $responder->($res->finalize);
                }
            );
            $pm->wait_all_children;
        };
    }
    
    method psgi_continuation_response (CodeRef $sub, HashRef $env) {
        my $page = $env->{PATH_INFO};
        
        my $app;
        unless ($self->continuation_page_stored($page)) {
            my $adapt = Continuity::Adapt::PSGI->new;
            my $continuity = Continuity->new(port => $self->port, adapter => $adapt, callback => $sub);
            $self->store_continuation($page => $adapt->loop_hook);
        }
        
        return &{ $self->continuation_for_page($page) }($env);
    }
    
    method transform_xml (Str $xml, Str $format) {
        my $source = $parser->load_xml(string => '<?xml version="1.0" encoding="ISO-8859-1"?><interface>' . $xml . '</interface>');
        my $stylesheet = $display_format_to_stylesheet{$format} || $display_format_to_stylesheet{html};
        my $result = $stylesheet->transform($source);
        return $stylesheet->output_as_chars($result);
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
    
    method log (Str $msg!, ArrayRef[Str] :$email_to?, Bool :$email_admin?, Str :$subject?, Str :$long_msg?, Bool :$force_when_testing?) {
        chomp($msg);
        
        # we could just warn to log to the log file if one is in use, but we'll
        # use flock to write to it for better multi-process behaviour
        my $log_msg = "$time{'yyyy/mm/dd hh:mm:ss'}: $msg\n";
        if ($self->_log_file_in_use) {
            my $log_file = $self->log_file;
            my $ok = open(my $fh, ">>", $log_file);
            if ($ok) {
                my $flock_fail = sub {
                    # failure to lock tends to be permanent, with the only
                    # apparent solution being to delete the log file; let's just
                    # assume the user is doing their own log rotation and so
                    # we're not losing too much stuff by this blind deletion...
                    close($fh);
                    unlink($log_file);
                    $ok = 0;
                };
                local $SIG{ALRM} = $flock_fail;
                alarm 60;
                $ok = flock($fh, LOCK_EX);
                alarm 0;
                if ($ok) {
                    seek($fh, 0, SEEK_END);
                    $ok = print $fh $log_msg;
                    flock($fh, LOCK_UN);
                }
                else {
                    &$flock_fail;
                }
                close($fh);
            }
            unless ($ok) {
                warn $log_msg, "Additionally, was unable to log to $log_file\n";
            }
        }
        else {
            warn $log_msg;
        }
        
        if (($force_when_testing || $self->deployment eq 'production') && ($email_to || $email_admin)) {
            # email the desired users
            my $domain      = $self->email_domain;
            my $admin_email = $self->admin_user . '@' . $domain;
            
            my $email = Email::Simple->create(
                header => [
                    To => $email_to ? join(', ', map { "$_\@$domain" } @$email_to) : $admin_email,
                    $email_admin && $email_to ? (Cc => $admin_email) : (),
                    From    => '"VRPipe Server" <vrpipe@do.not.reply>',
                    Subject => $subject || "VRPipe Server message",
                ],
                body => $msg . "\n" . ($long_msg || ''),
            );
            my $sent = Email::Sender::Simple->try_to_send($email);
            
            unless ($sent) {
                warn "$time{'yyyy/mm/dd hh:mm:ss'}: previous message failed to get sent to [", join(', ', ($email_to ? @$email_to : (), $email_admin ? '' : ())), "]\n";
            }
        }
    }
    
    method xml_tag (Str $tag, Str $cdata, Str $attribs?) {
        $attribs ||= '';
        $attribs &&= ' ' . $attribs;
        return '<' . $tag . $attribs . '><![CDATA[' . $cdata . ']]></' . $tag . '>';
    }
    
    method hash_to_xml (HashRef $hash, ArrayRef[Str] $key_order?) {
        $key_order ||= [sort { $a cmp $b } keys %$hash];
        
        my $xml = '<hash>';
        foreach my $key (@$key_order) {
            next unless defined $hash->{$key};
            $xml .= $self->xml_tag('pair', $hash->{$key}, qq[key="$key"]);
        }
        $xml .= '</hash>';
        
        return $xml;
    }
    
    method get_pipelinesetups (HashRef $opts, Bool $inactive?) {
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
        
        if ($multi_setups && !@setups) {
            die "No PipelineSetups match your settings (did you remember to specifiy --user?)\n";
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
}

1;
