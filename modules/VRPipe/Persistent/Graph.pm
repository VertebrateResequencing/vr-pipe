
=head1 NAME

VRPipe::Persistent::Graph - interface to a graph database

=head1 SYNOPSIS
    
    use VRPipe::Persistent::Graph;
    
    my $graph = VRPipe::Persistent::Graph->new();
    
    $graph->add_schema(
        namespace => 'QCGrind',
        label => 'Sample',    
        unique => [qw(sanger_id uuid)],
        indexed => [qw(public_name)]
    );
    
    my $node = $graph->add_node(
        namespace => 'QCGrind',
        label => 'Sample',
        properties => {
            sanger_id => 'sanger1',
            uuid => 'uuuuu',   
            public_name => 'public1'
        }
    );
    
    $graph->relate($node, $other_node, 'has');
    
    ($node) = $graph->get_nodes(
        namespace => 'QCGrind',
        label => 'Sample',
        properties => { public_name => 'public1' }
    );
    
    my ($related_node) = $graph->related_nodes(
        $node,
        namespace => 'QCGrind',
        label => 'Lane',
        max_depth => 4
    );

=head1 DESCRIPTION

For schema-less store of connected data we use a graph database; Neo4J in this
case.

This is essentially a wrapper around some cypher queries submitted to Neo4J via
its REST API, providing functions that can be used to store and retrieve
information about things.

Things (must) have a namespace, label and properties. A dynamically-applied
"schema" must be in place first, providing uniqueness constraints and indexes.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Persistent::Graph {
    use VRPipe::Config;
    use VRPipe::Persistent::SchemaBase;
    use Mojo::UserAgent;
    use JSON::XS;
    use Data::UUID;
    
    our $json       = JSON::XS->new->allow_nonref(1);
    our $data_uuid  = Data::UUID->new();
    our $vrp_config = VRPipe::Config->new();
    our ($ua, $transaction_endpoint, $global_label, $schemas, $schema_labels);
    our $ua_headers = { 'Accept' => 'application/json', 'Content-Type' => 'application/json', 'Charset' => 'UTF-8', 'X-Stream' => 'true' };
    
    has 'throw_with_no_stacktrace' => (
        is      => 'rw',
        isa     => 'Bool',
        default => 0
    );
    
    sub BUILD {
        my $self = shift;
        
        unless ($ua) {
            # we use Mojo::UserAgent instead of LWP::UserAgent because LWP has
            # some kind of truncation bug when we try to get very large
            # responses from Neo4J
            $ua = Mojo::UserAgent->new();
            
            # connect and get the transaction endpoint
            my $url = $vrp_config->neo4j_server_url();
            my $tx  = $ua->get($url => $ua_headers);
            my $res = $tx->success;
            unless ($res) {
                my $err = $tx->error;
                $self->throw("Failed to connect to '$url': [$err->{code}] $err->{message}");
            }
            my $decode = $json->decode($res->body);
            my $data_endpoint = $decode->{data} || $self->throw("No data endpoint found at $url");
            
            $tx = $ua->get($data_endpoint => $ua_headers);
            $res = $tx->success;
            unless ($res) {
                my $err = $tx->error;
                $self->throw("Failed to connect to '$data_endpoint': [$err->{code}] $err->{message}");
            }
            $decode = $json->decode($res->body);
            $transaction_endpoint = $decode->{transaction} || $self->throw("No transaction endpoint found at $data_endpoint");
            $transaction_endpoint .= '/commit';
            
            my $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
            if ($deployment eq 'production') {
                $global_label = "vdp";
            }
            else {
                my $user = getlogin || getpwuid($<);
                $global_label = "vdt$user";
            }
            $schema_labels = qq[`$global_label`:`Schema`];
        }
    }
    
    around throw (Str $msg) {
        if ($self->throw_with_no_stacktrace) {
            die $msg, "\n";
        }
        else {
            $self->$orig($msg);
        }
    }
    
    sub _run_cypher {
        my ($self, $array, $args) = @_;
        my $return_schema_nodes = $args->{return_schema_nodes};
        
        my $post_content = { statements => [] };
        foreach (@$array) {
            my ($cypher, $params) = @$_;
            push(
                @{ $post_content->{statements} },
                {
                    statement => $cypher,
                    $params ? (parameters => $params) : (),
                    resultDataContents => ['graph']
                }
            );
        }
        
        my $tx = $ua->post($transaction_endpoint => $ua_headers => json => $post_content);
        my $res = $tx->success;
        unless ($res) {
            my $err = $tx->error;
            $self->throw('[' . $err->{code} . '] ' . $err->{message});
        }
        my $decode = $json->decode($res->body);
        
        my $errors = $decode->{errors};
        if (@$errors) {
            my $error = $errors->[0];
            $self->throw('[' . $error->{code} . '] ' . $error->{message});
            #*** should auto-retry when the code matches TransientError
        }
        
        my (%nodes, %relationships);
        my $label_regex = qr/^$global_label\|([^\|]+)\|(.+)/;
        foreach my $result (@{ $decode->{results} }) {
            my $data = $result->{data} || next;
            foreach my $hash (@$data) {
                my $graph = $hash->{graph} || next;
                foreach my $node_details (@{ $graph->{nodes} || [] }) {
                    my $node_id = $node_details->{id};
                    next if exists $nodes{$node_id};
                    
                    # for speed reasons we don't create objects for nodes but just
                    # return the a hash and provide methods to extract stuff from
                    # the hash
                    
                    # convert labels to namespace and label
                    my ($namespace, $label);
                    my $ok = 0;
                    foreach my $this_label (@{ $node_details->{labels} }) {
                        if ($this_label =~ /$label_regex/) {
                            $namespace = $1;
                            $label     = $2;
                            $ok        = 1;
                            last;
                        }
                        elsif ($return_schema_nodes && $this_label eq $global_label) {
                            # Schema nodes
                            $ok = 1;
                        }
                    }
                    $ok || next; # only return nodes created by us
                    
                    my $node = { id => $node_id, properties => $node_details->{properties}, namespace => $namespace, label => $label };
                    $nodes{$node_id} = $node;
                }
                
                foreach my $rel_details (@{ $graph->{relationships} || [] }) {
                    my $rel_id = $rel_details->{id};
                    next if exists $relationships{$rel_id};
                    
                    # skip relationships if we skipped one of its nodes
                    next unless (exists $nodes{ $rel_details->{startNode} } && exists $nodes{ $rel_details->{endNode} });
                    
                    # again for speed reasons we just return the raw hash; this is
                    # more useful that applying the details to the nodes since this
                    # is the format needed for graph display
                    $relationships{$rel_id} = $rel_details;
                }
            }
        }
        
        if (defined wantarray()) {
            return { nodes => [values %nodes], relationships => [values %relationships] };
        }
        return;
    }
    
    method drop_database {
        $self->throw("drop_database() can only be used when testing") unless $global_label =~ /^vdt/;
        
        # drop all schemas (which drops all constraints and indexes)
        my @schema_nodes = @{ $self->_run_cypher([["MATCH (n:$schema_labels) RETURN n"]])->{nodes} };
        foreach my $node (@schema_nodes) {
            my $schema = $self->node_property($node, 'schema');
            my (undef, $namespace, $label) = split(/\|/, $schema);
            $self->drop_schema(namespace => $namespace, label => $label);
        }
        
        # drop all nodes and relationships
        $self->_run_cypher([["MATCH (n:`$global_label`) OPTIONAL MATCH (n:`$global_label`)-[r]-() DELETE n,r"]]);
        
        return 1;
    }
    
    sub _deployment_specific_label {
        my ($self, $namespace, $label) = @_;
        return "$global_label|$namespace|$label";
    }
    
    method add_schema (Str :$namespace!, Str :$label!, ArrayRef[Str] :$unique!, ArrayRef[Str] :$indexed?, ArrayRef[Str] :$required?) {
        # namespace and label cannot contain |
        foreach ($namespace, $label) {
            if (index($_, '|') != -1) {
                $self->throw("neither namespace or label may contain the | character");
            }
        }
        
        # have we already done this?
        my $dsl = $self->_deployment_specific_label($namespace, $label);
        my ($done) = @{ $self->_run_cypher([["MATCH (n:$schema_labels { schema: '$dsl' }) RETURN n"]], { return_schema_nodes => 1 })->{nodes} };
        unless ($done) {
            my @to_run;
            
            # set constraints (which also adds an index on the constraint)
            foreach my $field (@$unique) {
                if (index($field, '|') != -1) {
                    $self->throw("parameter may not contain the | character");
                }
                push(@to_run, ["CREATE CONSTRAINT ON (n:`$dsl`) ASSERT n.$field IS UNIQUE"]);
            }
            
            # add indexes
            foreach my $field (@{ $indexed || [] }) {
                if (index($field, '|') != -1) {
                    $self->throw("parameter may not contain the | character");
                }
                push(@to_run, ["CREATE INDEX ON :`$dsl`($field)"]);
            }
            
            $self->_run_cypher(\@to_run);
            # (sadly we can't do a single transaction that does both schema
            # updates above and the node creation below)
            
            # record that we've done this
            my $unique_fields = join('|', @$unique);
            my $indexed_arg = $indexed ? q[, indexed: '] . join('|', @$indexed) . q['] : '';
            my $required_arg = $required ? q[, required: '] . join('|', @$required) . q['] : '';
            $self->_run_cypher([["CREATE (:$schema_labels { schema: '$dsl', unique: '$unique_fields'$indexed_arg$required_arg })"]]);
            
            $schemas->{$dsl} = [$unique, $indexed || [], $required || []];
            
            return 1;
        }
        return 0;
    }
    
    method get_schema (Str :$namespace!, Str :$label!) {
        my $dsl = $self->_deployment_specific_label($namespace, $label);
        if (exists $schemas->{$dsl}) {
            return @{ $schemas->{$dsl} };
        }
        else {
            my ($schema) = @{ $self->_run_cypher([["MATCH (n:$schema_labels { schema: '$dsl' }) RETURN n"]], { return_schema_nodes => 1 })->{nodes} };
            if ($schema) {
                my $uniques  = [split(/\|/, $self->node_property($schema, 'unique'))];
                my $indexed  = [split(/\|/, $self->node_property($schema, 'indexed') || '')];
                my $required = [split(/\|/, $self->node_property($schema, 'required') || '')];
                $schemas->{$dsl} = [$uniques, $indexed, $required];
                return ($uniques, $indexed, $required);
            }
        }
    }
    
    method drop_schema (Str :$namespace!, Str :$label!) {
        my ($uniques, $indexed) = $self->get_schema(namespace => $namespace, label => $label);
        my $dsl = $self->_deployment_specific_label($namespace, $label);
        
        my @to_run;
        
        # remove constraints
        foreach my $field (@$uniques) {
            push(@to_run, ["DROP CONSTRAINT ON (n:`$dsl`) ASSERT n.$field IS UNIQUE"]);
        }
        
        # remove indexes
        foreach my $field (@$indexed) {
            push(@to_run, ["DROP INDEX ON :`$dsl`($field)"]);
        }
        
        # remove the node storing schema details, and our cache
        push(@to_run, ["MATCH (n:$schema_labels { schema: '$dsl' })-[r]-() DELETE n, r"]);
        
        $self->_run_cypher(\@to_run);
        
        delete $schemas->{$dsl};
    }
    
    sub _labels_and_param_map {
        my ($self, $namespace, $label, $params, $param_key, $check_required) = @_;
        
        # check that we have a schema for this
        my ($uniques, $indexed, $required) = $self->get_schema(namespace => $namespace, label => $label);
        $self->throw("You must first create a schema for namespace `$namespace` and label `$label`") unless $uniques;
        if ($check_required) {
            $self->throw("Parameters must be supplied") unless $params;
            foreach my $param (@$uniques, @{ $required || [] }) {
                $self->throw("Parameter '$param' must be supplied") unless defined $params->{$param};
            }
        }
        
        my $labels = "`$global_label`:`$global_label|$namespace|$label`";
        my $param_map = $params ? ' { ' . join(', ', map { "$_: {$param_key}.$_" } sort keys %$params) . ' }' : '';
        return ($labels, $param_map);
    }
    
    # in/outgoing HashRef is { type => 'type', node => $node }
    method add_nodes (Str :$namespace!, Str :$label!, ArrayRef[HashRef[Str]] :$properties!, HashRef :$incoming?, HashRef :$outgoing?) {
        my ($labels, $param_map) = $self->_labels_and_param_map($namespace, $label, $properties->[0], 'param', 1);
        
        my $left  = '';
        my $match = '';
        if ($incoming) {
            my $node_id = $self->node_id($incoming->{node});
            $left  = "(l)-[:$incoming->{type}]->";
            $match = "MATCH (l) WHERE id(l) = $node_id ";
        }
        my $right = '';
        if ($outgoing) {
            my $node_id = $self->node_id($outgoing->{node});
            $right = "-[:$outgoing->{type}]->(r)";
            $match .= "MATCH (r) WHERE id(r) = $node_id ";
        }
        
        my $cypher = "${match}MERGE $left(n:$labels$param_map)$right";
        $cypher .= ' RETURN n' if defined wantarray();
        
        my @to_run;
        foreach my $prop (@$properties) {
            push(@to_run, [$cypher, { 'param' => $prop }]);
        }
        
        if (defined wantarray()) {
            return @{ $self->_run_cypher(\@to_run)->{nodes} };
        }
        else {
            $self->_run_cypher(\@to_run);
        }
    }
    
    method add_node (Str :$namespace!, Str :$label!, HashRef[Str] :$properties!, HashRef :$incoming?, HashRef :$outgoing?) {
        my ($node) = $self->add_nodes(namespace => $namespace, label => $label, properties => [$properties], $incoming ? (incoming => $incoming) : (), $outgoing ? (outgoing => $outgoing) : ());
        return $node;
    }
    
    method delete_node (HashRef $node!) {
        $self->_run_cypher([["MATCH (n) WHERE id(n) = $node->{id} OPTIONAL MATCH (n)-[r]-() DELETE n, r"]]);
        return 1;
    }
    
    sub create_uuid {
        return $data_uuid->create_str();
    }
    
    method get_nodes (Str :$namespace!, Str :$label!, HashRef :$properties?) {
        my ($labels, $param_map) = $self->_labels_and_param_map($namespace, $label, $properties, 'param');
        return @{ $self->_run_cypher([["MATCH (n:$labels$param_map) RETURN n", { 'param' => $properties }]])->{nodes} };
    }
    
    method node_id (HashRef $node!) {
        if (defined $node->{id}) {
            return $node->{id};
        }
    }
    
    method node_namespace_and_label (HashRef $node!) {
        if (defined $node->{namespace} && defined $node->{label}) {
            return ($node->{namespace}, $node->{label});
        }
    }
    
    method node_property (HashRef $node!, Str $property!) {
        if (exists $node->{properties} && defined $node->{properties}->{$property}) {
            return $node->{properties}->{$property};
        }
    }
    
    method relate (HashRef $start_node!, HashRef $end_node!, Str :$type!) {
        return @{ $self->_run_cypher([["MATCH (a) WHERE id(a) = $start_node->{id} MATCH (b) WHERE id(b) = $end_node->{id} CREATE (a)-[r:$type]->(b) RETURN r"]])->{relationships} };
    }
    
    # incoming/outgoing/undirected hash refs are {min_depth, max_depth, type,
    # namespace, label, properties}, where the later 3 are result node specs and
    # with depths defaulting to 1 and others defaulting to undef; none supplied
    # defaults to undirected {min_depth => 1, max_depth => 1}
    # This returns hash of nodes and relationships for use by frontends;
    # related_nodes() calls this and returns just a list of nodes
    sub related {
        my ($self, $start_node, $undirected, $incoming, $outgoing, $result_nodes_only) = @_;
        $self->throw("undirected is mutually exclusive of in/outgoing") if $undirected && ($outgoing || $incoming);
        if (!$outgoing && !$incoming && !$undirected) {
            $undirected = { min_depth => 1, max_depth => 1 };
        }
        
        my $start_id = $self->node_id($start_node);
        if ($undirected) {
            my ($result_node_spec, $properties, $type, $min_depth, $max_depth) = $self->_related_nodes_hashref_parse($undirected, 'param');
            my $return = $result_nodes_only ? 'u' : 'p';
            return $self->_run_cypher([["MATCH p = (start)-[$type*$min_depth..$max_depth]-(u$result_node_spec) WHERE id(start) = $start_id RETURN $return", { 'param' => $properties }]]);
        }
        else {
            my (%all_properties, @return);
            my $left = '';
            if ($incoming) {
                my ($result_node_spec, $properties, $type, $min_depth, $max_depth) = $self->_related_nodes_hashref_parse($incoming, 'left');
                $left = "(l$result_node_spec)-[$type*$min_depth..$max_depth]->";
                push(@return, 'l');
                $all_properties{left} = $properties if $properties;
            }
            my $right = '';
            if ($outgoing) {
                my ($result_node_spec, $properties, $type, $min_depth, $max_depth) = $self->_related_nodes_hashref_parse($outgoing, 'right');
                $right = "-[$type*$min_depth..$max_depth]->(r$result_node_spec)";
                push(@return, 'r');
                $all_properties{right} = $properties if $properties;
            }
            
            my $return;
            if ($result_nodes_only) {
                $return = join(', ', @return);
            }
            else {
                $return = 'p';
            }
            return $self->_run_cypher([["MATCH p = $left(start)$right where id(start) = $start_id RETURN $return", keys %all_properties ? \%all_properties : ()]]);
        }
    }
    
    sub _related_nodes_hashref_parse {
        my ($self, $hashref, $param_key) = @_;
        my $type = $hashref->{type} ? ":`$hashref->{type}`" : '';
        my $min_depth = $hashref->{min_depth} || 1;
        my $max_depth = $hashref->{max_depth} || $min_depth;
        my $result_node_spec = '';
        if ($hashref->{namespace} && $hashref->{label}) {
            my ($labels, $param_map) = $self->_labels_and_param_map($hashref->{namespace}, $hashref->{label}, $hashref->{properties}, $param_key);
            $result_node_spec = ':' . $labels . $param_map;
        }
        return ($result_node_spec, $hashref->{properties}, $type, $min_depth, $max_depth);
    }
    
    method related_nodes (HashRef $start_node!, HashRef :$outgoing?, HashRef :$incoming?, HashRef :$undirected?) {
        return @{ $self->related($start_node, $undirected, $incoming, $outgoing, 1)->{nodes} };
    }
    
    method root_nodes {
        return @{ $self->_run_cypher([["MATCH (n:`$global_label`) OPTIONAL MATCH (n)<-[r]-() WITH n,r WHERE r IS NULL RETURN n"]])->{nodes} };
    }
    
    sub json_encode {
        shift;
        return $json->encode(shift);
    }
    
    sub json_decode {
        shift;
        return $json->decode(shift);
    }
}

1;
