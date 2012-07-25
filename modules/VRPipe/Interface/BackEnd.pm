
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
    
    my $xsl_html = <<XSL;
<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template match="/">
<html>
    <body>
        <xsl:apply-templates/>  
    </body>
</html>
</xsl:template>

<xsl:template match="title">
    <h2><xsl:value-of select="."/></h2>
</xsl:template>

<xsl:template match="table">
    <table border="1">
        <xsl:apply-templates select="headings"/>  
        <xsl:apply-templates select="row"/>
    </table>
</xsl:template>

<xsl:template match="headings">
    <tr bgcolor="#9acd32">
        <xsl:apply-templates select="heading"/>  
    </tr>
</xsl:template>

<xsl:template match="heading">
    <th><xsl:value-of select="."/></th>
</xsl:template>

<xsl:template match="row">
    <tr>
        <xsl:choose>
            <xsl:when test="type = 'PipelineSetup'">
                <xsl:choose>
                    <xsl:when test="active = 0">
                        <td style="color:#cccccc"><xsl:value-of select="id"/></td>
                        <td style="color:#cccccc"><xsl:value-of select="name"/></td>
                        <td style="color:#cccccc"><xsl:value-of select="user"/></td>
                        <td style="color:#cccccc"><xsl:value-of select="active"/></td>
                    </xsl:when>
                    <xsl:otherwise>
                        <td><xsl:value-of select="id"/></td>
                        <td><xsl:value-of select="name"/></td>
                        <td><xsl:value-of select="user"/></td>
                        <td><xsl:value-of select="active"/></td>
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:when>
            <xsl:otherwise>
                <td>unsuported type</td>
            </xsl:otherwise>
        </xsl:choose>
    </tr>
</xsl:template>
</xsl:stylesheet>
XSL
    my $xsl_plain = <<XSL;
<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text"/>

<xsl:template match="/">
    <xsl:apply-templates/>
</xsl:template>

<xsl:template match="table">
    <xsl:apply-templates select="headings"/>  
    <xsl:apply-templates select="row"/>
    <xsl:text>\n</xsl:text>
</xsl:template>

<xsl:template match="headings">
    <xsl:apply-templates select="heading"/><xsl:text>\n</xsl:text>
</xsl:template>

<xsl:template match="heading">
    <xsl:value-of select="."/><xsl:text>\t</xsl:text>
</xsl:template>

<xsl:template match="row">
    <xsl:choose>
        <xsl:when test="type = 'PipelineSetup'">
            <xsl:value-of select="id"/><xsl:text>\t</xsl:text>
            <xsl:value-of select="name"/><xsl:text>\t</xsl:text>
            <xsl:value-of select="user"/><xsl:text>\t</xsl:text>
            <xsl:value-of select="active"/><xsl:text>\t\n</xsl:text>
        </xsl:when>
        <xsl:otherwise>
            <xsl:text>unsuported type\n</xsl:text>
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>

</xsl:stylesheet>
XSL
    our %display_format_to_xsl = (html  => $xsl_html,
                                  plain => $xsl_plain);
    
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
        
        unless (defined $opts{display_format}) {
            $opts{display_format} = 'html';
        }
        
        return \%opts;
    }
    
    method output ($req, $xml, $format) {
        my $parser = XML::LibXML->new();
        my $xslt   = XML::LibXSLT->new();
        
        my $source = $parser->load_xml(string => '<?xml version="1.0" encoding="ISO-8859-1"?>' . $xml);
        my $style_doc = $parser->load_xml(string => $display_format_to_xsl{$format} || $self->throw("invalid display_format '$format'"));
        
        my $stylesheet = $xslt->parse_stylesheet($style_doc);
        my $results    = $stylesheet->transform($source);
        
        $req->respond({ content => ["text/$format", $stylesheet->output_string($results)] });
    }
    
    method display_hash (Str $name, HashRef $hash, ArrayRef[Str] $key_order?) {
        $key_order ||= [sort { $a cmp $b } keys %$hash];
        my $xml = "$name:";
        my ($extra_tabs) = $name =~ /^(\t+)/;
        $extra_tabs ||= '';
        foreach my $key (@$key_order) {
            next unless defined $hash->{$key};
            $xml .= "$extra_tabs\t" . $key . ' => ' . $hash->{$key};
        }
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
