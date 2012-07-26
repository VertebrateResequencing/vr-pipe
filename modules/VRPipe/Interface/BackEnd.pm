
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
    use AnyEvent::HTTPD;
    use Sys::Hostname::Long;
    use VRPipe::Config;
    use XML::LibXSLT;
    use XML::LibXML;
    
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

<xsl:template match="/interface/title">
    <h1><xsl:value-of select="."/></h1>
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

<xsl:template match="/interface/title"/>

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
    
    <xsl:if test="@display_mode!='list'">
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
    our %display_format_to_stylesheet = (html  => $xslt->parse_stylesheet($parser->load_xml(string => $xsl_html)),
                                         plain => $xslt->parse_stylesheet($parser->load_xml(string => $xsl_plain)));
    
    has 'port' => (is       => 'ro',
                   isa      => PositiveInt,
                   required => 1);
    
    has 'schema' => (is      => 'ro',
                     isa     => 'VRPipe::Persistent::Schema',
                     lazy    => 1,
                     builder => '_build_schema');
    
    has 'httpd' => (is      => 'ro',
                    isa     => 'AnyEvent::HTTPD',
                    lazy    => 1,
                    builder => '_build_httpd');
    
    method _build_schema {
        my $m = VRPipe::Manager->get;
        return $m->result_source->schema;
    }
    
    method _build_httpd {
        my $port = $self->port;
        warn "The web interface can be reached at http://", hostname_long, ":$port/\n";
        return AnyEvent::HTTPD->new(port => $port); # host `uname -n`; parse the ip address, use as host? Currently defaults to 0.0.0.0
    }
    
    sub BUILD {
        my $self = shift;
        
        # the idea is (per host) we'll have 1 daemon with 1 BackEnd instance per
        # dsn, and each dsn is supposed to have its own port associated with it.
        # we look at the port and figure out if it corresponds to the current
        # user's production or testing database; if neither, we die
        my $port       = $self->port;
        my $vrp_config = VRPipe::Config->new();
        my $deployment;
        foreach my $dep (qw(production testing)) {
            my $method_name = $dep . '_interface_port';
            my $this_port   = $vrp_config->$method_name();
            if ($this_port == $port) {
                $deployment = $dep;
            }
        }
        unless ($deployment) {
            $self->throw("VRPipe interface port $port did not correspond to a deployment");
        }
        
        VRPipe::Persistent::SchemaBase->database_deployment($deployment);
        require VRPipe::Persistent::Schema;
    }
    
    method req_to_opts ($req, ArrayRef $ctp_array?) {
        my %opts = $req->vars;
        
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
                            $self->die_with_error("$desired is not a valid $class id");
                        }
                    }
                    else {
                        ($found) = $full_class->search({ name => $desired });
                        unless ($found) {
                            $self->die_with_error("$desired is not a valid $class name");
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
                        $self->error("--$opt @desired resulted in more than one $class object; only using the first");
                    }
                }
            }
        }
        
        return \%opts;
    }
    
    method output ($req, $xml) {
        my $format = $req->parm('display_format');
        $format ||= 'html';
        
        my $source = $parser->load_xml(string => '<?xml version="1.0" encoding="ISO-8859-1"?>' . $xml);
        my $stylesheet = $display_format_to_stylesheet{$format} || $self->throw("invalid display_format '$format'");
        my $result = $stylesheet->transform($source);
        
        $req->respond({ content => ["text/$format", $stylesheet->output_as_chars($result)] });
    }
    
    method hash_to_xml (HashRef $hash, ArrayRef[Str] $key_order?) {
        $key_order ||= [sort { $a cmp $b } keys %$hash];
        
        my $xml = '<hash>';
        foreach my $key (@$key_order) {
            next unless defined $hash->{$key};
            $xml .= qq[<pair key="$key"><![CDATA[$hash->{$key}]]></pair>];
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
            $self->die_with_error("No PipelineSetups match your settings (did you remember to specifiy --user?)");
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
