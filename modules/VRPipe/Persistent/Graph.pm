
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

Copyright (c) 2014, 2015 Genome Research Limited.

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
    use VRPipe::Persistent::InMemory;
    use Mojo::UserAgent;
    use JSON::XS;
    use Data::UUID;
    use DateTime::Format::Natural;
    use MIME::Base64;
    use URI::Escape;
    
    our $json       = JSON::XS->new->allow_nonref(1);
    our $data_uuid  = Data::UUID->new();
    our $vrp_config = VRPipe::Config->new();
    our ($ua, $url, $transaction_endpoint, $global_label, $schemas, $schema_labels);
    our $ua_headers = { 'Accept' => 'application/json', 'Content-Type' => 'application/json', 'Charset' => 'UTF-8', 'X-Stream' => 'true' };
    
    has 'throw_with_no_stacktrace' => (
        is      => 'rw',
        isa     => 'Bool',
        default => 0
    );
    
    has '_collected' => (
        is      => 'ro',
        traits  => ['Array'],
        isa     => 'ArrayRef',
        default => sub { [] },
        handles => {
            '_collect'          => 'push',
            '_empty_collection' => 'clear'
        }
    );
    
    sub BUILD {
        my $self = shift;
        
        unless ($ua) {
            my $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
            
            # we use Mojo::UserAgent instead of LWP::UserAgent because LWP has
            # some kind of truncation bug when we try to get very large
            # responses from Neo4J
            $ENV{MOJO_MAX_MESSAGE_SIZE} = 0; # avoid Maximum message size exceeded errors when we get lots of data from a query
            $ua = Mojo::UserAgent->new();
            $ua->connect_timeout(60)->inactivity_timeout(0)->request_timeout(0);
            
            # connect and get the transaction endpoint
            my $method_name = $deployment . '_neo4j_server_url';
            $url         = $vrp_config->$method_name();
            $method_name = $deployment . '_neo4j_user';
            my $user = $vrp_config->$method_name();
            $method_name = $deployment . '_neo4j_password';
            my $password = $vrp_config->$method_name();
            
            if ($user && $password) {
                # Requests should include an Authorization header, with a value
                # of Basic <payload>, where "payload" is a base64 encoded string
                # of "username:password"
                $ua_headers->{Authorization} = 'Basic ' . substr(encode_base64("$user:$password"), 0, -2);
            }
            
            foreach my $try_num (1 .. 20) {
                eval {
                    my $tx = $ua->get("$url" => $ua_headers);
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
                };
                if ($@) {
                    if ($try_num == 20) {
                        die $@;
                    }
                    else {
                        sleep($try_num * 2);
                    }
                }
                else {
                    last;
                }
            }
            
            if ($deployment eq 'production') {
                $global_label = "vdp";
            }
            else {
                my $user = getpwuid($<);
                $user || $self->throw("Could not determine user, so can't write anything to the graph database");
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
    
    sub _global_label {
        return $global_label;
    }
    
    sub _run_cypher {
        my ($self, $array, $args) = @_;
        my $return_schema_nodes  = $args->{return_schema_nodes};
        my $return_history_nodes = $args->{return_history_nodes};
        
        my $post_content = { statements => [] };
        my $example_cypher;
        foreach (@$array) {
            my ($cypher, $params) = @$_;
            $example_cypher = $cypher;
            push(
                @{ $post_content->{statements} },
                {
                    statement => $cypher,
                    $params ? (parameters => $params) : (),
                    resultDataContents => ['graph']
                }
            );
        }
        
        my $decode;
        foreach my $try_num (1 .. 20) {
            my ($tx, $res);
            eval { $tx = $ua->post($transaction_endpoint => $ua_headers => json => $post_content); };
            $res = $tx->success if $tx;
            unless ($res) {
                my $err  = $tx->error;
                my $code = $err->{code};
                $code ||= 'no error code';
                my $message = $err->{message};
                $message ||= $@ || '(no message)';
                
                if ($message =~ /connection/i) {
                    # neo4j server may be down during a backup, so we'll wait a
                    # a few mins for it to come back
                    warn "retrying cypher [$example_cypher] due to: [$code] $message\n";
                    sleep(2 * $try_num);
                    next;
                }
                
                $self->throw("[$code] [$example_cypher] $message");
            }
            $decode = $json->decode($res->body);
            
            my $errors = $decode->{errors};
            if (@$errors) {
                my $err  = $errors->[0];
                my $code = $err->{code};
                $code ||= 'no error code';
                my $message = $err->{message};
                $message ||= '(no message)';
                if ($try_num < 20 && ($code =~ /CouldNotCommit|TransientError/ || $message =~ /connection/i)) {
                    warn "retrying cypher [$example_cypher] due to: [$code] $message\n";
                    
                    if ($message =~ /connection/i) {
                        sleep(2 * $try_num);
                    }
                    else {
                        sleep(1);
                    }
                }
                else {
                    $self->throw("[$code] [$example_cypher] $message");
                }
            }
            else {
                last;
            }
        }
        
        my (%nodes, $node_order, %relationships);
        my $label_regex = qr/^$global_label\|([^\|]+)\|(.+)/;
        foreach my $result (@{ $decode->{results} }) {
            my $data = $result->{data} || next;
            foreach my $hash (@$data) {
                my $graph = $hash->{graph} || next;
                foreach my $node_details (@{ $graph->{nodes} || [] }) {
                    my $node_id = $node_details->{id};
                    next if exists $nodes{$node_id};
                    
                    # for speed reasons we don't create objects for nodes but
                    # just return the hash and provide methods to extract stuff
                    # from the hash
                    
                    # convert labels to namespace and label
                    my ($namespace, $label);
                    my $ok = 0;
                    foreach my $this_label (@{ $node_details->{labels} }) {
                        if ($this_label =~ /$label_regex/) {
                            $namespace = $1;
                            
                            unless ($return_history_nodes) {
                                if ($namespace eq 'PropertiesWithHistory') {
                                    last;
                                }
                            }
                            
                            $label = $2;
                            $ok    = 1;
                            last;
                        }
                        elsif ($return_schema_nodes && $this_label eq $global_label) {
                            # Schema nodes
                            $ok = 1;
                        }
                    }
                    $ok || next; # only return nodes created by us
                    
                    my $node = { id => int($node_id), properties => $node_details->{properties}, namespace => $namespace, label => $label };
                    $nodes{$node_id} = [++$node_order, $node];
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
            return { nodes => [map { $_->[1] } sort { $a->[0] <=> $b->[0] } values %nodes], relationships => [values %relationships] };
        }
        return;
    }
    
    method drop_database {
        $self->throw("drop_database() can only be used when testing") unless $global_label =~ /^vdt/;
        
        # drop all schemas (which drops all constraints and indexes)
        my @schema_nodes = @{ $self->_run_cypher([["MATCH (n:$schema_labels) RETURN n"]], { return_schema_nodes => 1 })->{nodes} };
        foreach my $node (@schema_nodes) {
            my $schema = $self->node_property($node, 'schema');
            my (undef, $namespace, $label) = split(/\|/, $schema);
            $self->drop_schema(namespace => $namespace, label => $label);
        }
        
        # drop all nodes and relationships
        $self->_run_cypher([["MATCH (n:`$global_label`) OPTIONAL MATCH (n)-[r]-() DELETE n,r"]]);
        
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
        
        # have we already done this? block until locked so multiple processes
        # don't end up corrupting a possible update
        my $im       = VRPipe::Persistent::InMemory->new();
        my $lock_key = 'graph.' . $namespace . '.' . $label . '.schema_update';
        $im->block_until_locked($lock_key);
        my ($current_unique, $current_indexed, $current_required) = $self->get_schema(namespace => $namespace, label => $label);
        
        if ($current_unique) {
            # see if it's changed and we need to update
            my $changed = 0;
            CMP: foreach my $cmp ([$unique, $current_unique], [$indexed, $current_indexed], [$required, $current_required]) {
                my ($new, $current) = @$cmp;
                if ($new && !@$new) {
                    undef $new;
                }
                if ($current && !@$current) {
                    undef $current;
                }
                
                if (($new && !$current) || (!$new && $current)) {
                    $changed = 1;
                    last;
                }
                next unless $new;
                
                my %current = map { $_ => 1 } @$current;
                foreach my $param (@$new) {
                    unless (exists $current{$param}) {
                        $changed = 1;
                        last CMP;
                    }
                    delete $current{$param};
                }
                
                if (keys %current) {
                    $changed = 1;
                    last;
                }
            }
            
            if ($changed) {
                $self->drop_schema(namespace => $namespace, label => $label);
                undef $current_unique;
            }
        }
        
        unless ($current_unique) {
            my @to_run;
            
            # set constraints (which also adds an index on the constraint)
            my $dsl = $self->_deployment_specific_label($namespace, $label);
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
            $self->_run_cypher([["MERGE (:$schema_labels { schema: '$dsl', unique: '$unique_fields'$indexed_arg$required_arg })"]]);
            
            $schemas->{$dsl} = [$unique, $indexed || [], $required || []];
            $im->unlock($lock_key);
            return 1;
        }
        
        $im->unlock($lock_key);
        return 0;
    }
    
    sub get_schema {
        my ($self, %opts) = @_;
        my $namespace = $opts{namespace};
        my $label     = $opts{label};
        
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
        
        # (these need to be run separately from the following)
        $self->_run_cypher(\@to_run);
        
        # remove the node storing schema details, and our cache
        $self->_run_cypher([["MATCH (n:$schema_labels { schema: '$dsl' }) OPTIONAL MATCH (n)-[r]-() DELETE n, r"]]);
        
        delete $schemas->{$dsl};
    }
    
    sub _labels {
        my ($self, $namespace, $label) = @_;
        return "`$global_label`:`$global_label|$namespace|$label`";
    }
    
    sub _param_map {
        my ($self, $params, $param_key) = @_;
        return $params ? ' { ' . join(', ', map { "`$_`: {$param_key}.`$_`" } sort keys %$params) . ' }' : '';
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
        
        return ($self->_labels($namespace, $label), $self->_param_map($params, $param_key));
    }
    
    # in/outgoing HashRef is { type => 'type', node => $node }
    # instead of node you can have node_spec => { namespace => ..., label => ..., properties => { ... } }
    # to have multiple in/outgoing, supply an array ref of those hashrefs
    sub add_nodes {
        my ($self, %opts) = @_;
        my $namespace            = $opts{namespace};
        my $label                = $opts{label};
        my $properties           = $opts{properties};
        my $update               = $opts{update} || 0;
        my $return_history_nodes = $opts{return_history_nodes} || 0;
        my $incoming             = $opts{incoming};
        my $outgoing             = $opts{outgoing};
        my $enqueue              = $opts{enqueue} || 0;
        
        my ($labels, $param_map);
        my $set = '';
        if ($update) {
            # split out unique params from the others; we'll merge on the
            # uniques and set the remainder
            my ($uniques, $indexed, $required) = $self->get_schema(namespace => $namespace, label => $label);
            my %uniques  = map { $_ => 1 } @$uniques;
            my %required = map { $_ => 1 } @$required;
            my ($unique_props, $other_props);
            while (my ($key, $val) = each %{ $properties->[0] }) {
                if (exists $uniques{$key}) {
                    $unique_props->{$key} = $val;
                }
                else {
                    $other_props->{$key} = $val;
                }
                
                delete $uniques{$key};
                delete $required{$key};
            }
            
            foreach my $param (keys %uniques, keys %required) {
                $self->throw("Parameter '$param' must be supplied");
            }
            
            if (keys %$unique_props && keys %$other_props) {
                ($labels, $param_map) = $self->_labels_and_param_map($namespace, $label, $unique_props, 'param', 0);
                $set = ' SET n += ' . $self->_param_map($other_props, 'param');
            }
        }
        unless ($labels) {
            ($labels, $param_map) = $self->_labels_and_param_map($namespace, $label, $properties->[0], 'param', 1);
        }
        
        my $rel_merge    = '';
        my $match        = '';
        my $extra_params = {};
        if ($incoming) {
            ($rel_merge, $match) = $self->_add_rel_cypher($incoming, 'l', 'in', $extra_params);
        }
        if ($outgoing) {
            my ($out_merge, $out_match) = $self->_add_rel_cypher($outgoing, 'r', 'out', $extra_params);
            $rel_merge .= $out_merge;
            $match     .= $out_match;
        }
        
        my $cypher = "${match}MERGE (n:$labels$param_map)$set$rel_merge";
        $cypher .= ' RETURN n' if defined wantarray();
        
        my @to_run;
        foreach my $prop (@$properties) {
            push(@to_run, [$cypher, { 'param' => $prop, %{$extra_params} }]);
        }
        
        if ($enqueue) {
            $self->_collect(@to_run);
            return;
        }
        
        if (defined wantarray()) {
            return @{ $self->_run_cypher(\@to_run, { return_history_nodes => $return_history_nodes })->{nodes} };
        }
        else {
            $self->_run_cypher(\@to_run);
        }
    }
    
    method dispatch_queue {
        my @nodes;
        if (defined wantarray()) {
            @nodes = @{ $self->_run_cypher($self->_collected, { return_history_nodes => 0 })->{nodes} };
        }
        else {
            $self->_run_cypher($self->_collected);
        }
        
        $self->_empty_collection;
        
        return @nodes;
    }
    
    sub _add_rel_cypher {
        my ($self, $spec, $prefix, $direction, $extra_params) = @_;
        my $i     = 0;
        my $merge = '';
        my $match = '';
        foreach my $s (ref($spec) eq 'ARRAY' ? @$spec : ($spec)) {
            $i++;
            
            my $type = $s->{type} || $self->throw("type must be defined in incoming/outgoing");
            
            if ($direction eq 'in') {
                $merge .= " MERGE (`$prefix$i`)-[:$type]->(n)";
            }
            elsif ($direction eq 'out') {
                $merge .= " MERGE (n)-[:$type]->(`$prefix$i`)";
            }
            else {
                $self->throw("bad direction '$direction' supplied to _add_rel_cypher");
            }
            
            if ($s->{node}) {
                my $node_id = $self->node_id($s->{node});
                $match .= "MATCH (`$prefix$i`) WHERE id(`$prefix$i`) = $node_id ";
            }
            elsif ($s->{node_spec}) {
                my $node_spec = $s->{node_spec};
                my $param_key = "$prefix${i}param";
                my ($labels, $param_map) = $self->_labels_and_param_map($node_spec->{namespace}, $node_spec->{label}, $node_spec->{properties}, $param_key);
                $match .= "MATCH (`$prefix$i`:$labels$param_map)";
                $extra_params->{$param_key} = $node_spec->{properties};
            }
            else {
                $self->throw("node or node_spec must be defined in incoming/outgoing");
            }
        }
        return ($merge, $match);
    }
    
    sub add_node {
        my ($self, %opts) = @_;
        my $namespace  = $opts{namespace};
        my $label      = $opts{label};
        my $properties = $opts{properties};
        my $update     = $opts{update} || 0;
        my $incoming   = $opts{incoming};
        my $outgoing   = $opts{outgoing};
        my $enqueue    = $opts{enqueue} || 0;
        
        my ($node) = $self->add_nodes(namespace => $namespace, label => $label, properties => [$properties], enqueue => $enqueue, update => $update, $incoming ? (incoming => $incoming) : (), $outgoing ? (outgoing => $outgoing) : ());
        return $node;
    }
    
    method delete_node (HashRef|Object $node!) {
        $self->_run_cypher([["MATCH (n) WHERE id(n) = $node->{id} OPTIONAL MATCH (n)-[r]-() DELETE n, r"]]);
        return 1;
    }
    
    sub create_uuid {
        return $data_uuid->create_str();
    }
    
    method get_nodes (Str :$namespace!, Str :$label!, Bool :$return_history_nodes = 0, HashRef :$properties?) {
        my ($labels, $param_map) = $self->_labels_and_param_map($namespace, $label, $properties, 'param');
        return @{ $self->_run_cypher([["MATCH (n:$labels$param_map) RETURN n", { 'param' => $properties }]], { return_history_nodes => $return_history_nodes })->{nodes} };
    }
    
    method get_node_by_id (Int $id) {
        my @nodes = @{ $self->_run_cypher([["MATCH (n) WHERE id(n) = $id RETURN n"]], { return_history_nodes => 1 })->{nodes} };
        return $nodes[0];
    }
    
    # properties is { property => 'regex', ... }
    method get_nodes_by_regex (Str :$namespace!, Str :$label!, HashRef :$properties!) {
        my $labels = $self->_labels($namespace, $label);
        
        my %regexes;
        my @wheres;
        while (my ($key, $val) = each %$properties) {
            $regexes{$key} = $val;
            my $operator = $val =~ /^[\w\-#]+$/ ? '=' : '=~';
            push(@wheres, "n.`$key` $operator {param}.`$key`");
        }
        my $where = join(' AND ', @wheres);
        
        return @{ $self->_run_cypher([[qq{MATCH (n:$labels) WHERE $where RETURN n}, { 'param' => \%regexes }]])->{nodes} };
    }
    
    sub node_id {
        my ($self, $node) = @_;
        if (defined $node->{id}) {
            return $node->{id};
        }
    }
    
    method node_namespace_and_label (HashRef|Object $node!) {
        if (defined $node->{namespace} && defined $node->{label}) {
            return ($node->{namespace}, $node->{label});
        }
    }
    
    method node_properties (HashRef|Object $node!, Bool :$flatten_parents = 0) {
        if ($flatten_parents) {
            unless (exists $node->{parent_properties}) {
                # get all the node properties of all parent nodes
                $node->{parent_properties} = {};
                foreach my $parent ($self->related_nodes($node, incoming => { min_depth => 1, max_depth => 6 })) { # we can't have max_depth 999 since that can kill Neo4J
                    my $prefix = $parent->{namespace} ne $node->{namespace} ? $parent->{namespace} . '_' : '';
                    $prefix .= $parent->{label};
                    $prefix = lc($prefix);
                    
                    while (my ($key, $val) = each %{ $parent->{properties} || {} }) {
                        my $full_key = $prefix . '_' . $key;
                        $node->{parent_properties}->{$full_key} = $val;                                            #*** we should probably merge and keep multiple vals for same full_key
                    }
                }
            }
            
            return { %{ $node->{parent_properties} }, %{ $node->{properties} || {} } };
        }
        
        return $node->{properties} || {};
    }
    
    method node_property (HashRef|Object $node!, Str $property!, Bool :$check_parents = 0) {
        my $properties = $self->node_properties($node, flatten_parents => $check_parents);
        return $properties->{$property} if exists $properties->{$property};
        return;
    }
    
    method node_add_properties (HashRef|Object $node!, HashRef $properties!) {
        my $id = $self->node_id($node);
        my $properties_map = $self->_param_map($properties, 'param');
        # (this requires Neo4J v2.1.2 +)
        my ($updated_node) = @{ $self->_run_cypher([["MATCH (n) WHERE id(n) = $id SET n += $properties_map return n", { 'param' => $properties }]])->{nodes} };
        $node->{properties} = $updated_node->{properties};
        return;
    }
    
    method node_set_properties (HashRef|Object $node!, HashRef $properties!) {
        my $id = $self->node_id($node);
        my ($updated_node) = @{ $self->_run_cypher([["MATCH (n) WHERE id(n) = $id SET n = { props } return n", { 'props' => $properties }]])->{nodes} };
        $node->{properties} = $updated_node->{properties};
        return;
    }
    
    method node_remove_property (HashRef|Object $node!, Str $property!) {
        my $id = $self->node_id($node);
        my ($updated_node) = @{ $self->_run_cypher([["MATCH (n) WHERE id(n) = $id REMOVE n.`$property` return n"]])->{nodes} };
        $node->{properties} = $updated_node->{properties};
        return;
    }
    
    # selfish => 1 means that the start node can be related to unlimited
    # end_nodes, but the end_node can only be related to a single node with the
    # same label (and relationship type) as the start node.
    # replace => 1 means that the end node can be related to unlimited
    # start_nodes, but the start_node can only be related to a single end_node
    # with the same label (and relationship type) as the end_node.
    # selfish => 1, replace => 1 gives you a 1:1 relationship between nodes of
    # the relevant labels and with the given relationship type.
    # properties option are properties to set on the relationship itself.
    method relate (HashRef|Object $start_node!, HashRef|Object $end_node!, Str :$type!, HashRef :$properties?, Bool :$selfish = 0, Bool :$replace = 0) {
        my @cypher;
        
        if ($selfish) {
            my $labels = $self->_labels($start_node->{namespace}, $start_node->{label});
            push(@cypher, ["MATCH (a)<-[r:$type]-(b:$labels) WHERE id(a) = $end_node->{id} AND id(b) <> $start_node->{id} DELETE r"]);
        }
        
        if ($replace) {
            my $labels = $self->_labels($end_node->{namespace}, $end_node->{label});
            push(@cypher, ["MATCH (a)-[r:$type]->(b:$labels) WHERE id(a) = $start_node->{id} AND id(b) <> $end_node->{id} DELETE r"]);
        }
        
        my $r_props = '';
        if ($properties) {
            my $map = $self->_param_map($properties, 'param');
            $r_props = " SET r = $map";
        }
        
        push(@cypher, ["MATCH (a),(b) WHERE id(a) = $start_node->{id} AND id(b) = $end_node->{id} MERGE (a)-[r:$type]->(b)$r_props RETURN r", $r_props ? ({ 'param' => $properties }) : ()]);
        
        return @{ $self->_run_cypher(\@cypher)->{relationships} };
    }
    
    # each hashref in $spec_list is { from => { id }|{ namespace, label, properties }, to => {...}, type => '...', properties => {} }
    method create_mass_relationships (ArrayRef[HashRef] $spec_list) {
        my @cypher;
        
        foreach my $spec (@$spec_list) {
            my $type = $spec->{type} || $self->throw("type must be present");
            
            my $match  = '';
            my $params = {};
            foreach my $direction ('from', 'to') {
                my $node_spec = $spec->{$direction} || $self->throw("$direction must be present");
                if ($node_spec->{id}) {
                    $match .= "MATCH ($direction) WHERE id($direction) = $node_spec->{id} ";
                }
                elsif ($node_spec->{namespace} && $node_spec->{label} && $node_spec->{properties}) {
                    my ($labels, $param_map) = $self->_labels_and_param_map($node_spec->{namespace}, $node_spec->{label}, $node_spec->{properties}, $direction);
                    $match .= "MATCH ($direction:$labels$param_map) ";
                    $params->{$direction} = $node_spec->{properties};
                }
                else {
                    $self->throw("node id or namespace,label,properties must be defined in $direction");
                }
            }
            
            my $r_props = '';
            if ($spec->{properties}) {
                my $map = $self->_param_map($spec->{properties}, 'relparams');
                $r_props = " SET r = $map";
                $params->{relparams} = $spec->{properties};
            }
            
            push(@cypher, [$match . " MERGE (from)-[r:$type]->(to)$r_props", $params]);
        }
        
        $self->_run_cypher(\@cypher);
    }
    
    # delete all relationships between 2 nodes, optionally limited by type
    method divorce (HashRef|Object $start_node!, HashRef|Object $end_node!, Str :$type?) {
        $self->mass_divorce([[$start_node, $end_node, $type ? ($type) : ()]]);
    }
    
    # each member of $spec_list is [$start_node, $end_node, $type], where $type
    # is optional
    method mass_divorce (ArrayRef $spec_list!) {
        my @cypher;
        
        foreach my $spec (@$spec_list) {
            my ($start_node, $end_node, $type) = @$spec;
            $type ||= '';
            $type &&= ':' . $type;
            push(@cypher, ["MATCH (a)-[rel$type]-(b) WHERE id(a) = $start_node->{id} AND id(b) = $end_node->{id} DELETE rel"]);
        }
        
        $self->_run_cypher(\@cypher);
    }
    
    method relationship_set_properties (HashRef $rel!, HashRef $properties!) {
        my $id = $rel->{id};
        my $properties_map = $self->_param_map($properties, 'param');
        my ($updated_rel) = @{ $self->_run_cypher([["MATCH ()-[r]->() WHERE id(r) = $id SET r = $properties_map return r", { 'param' => $properties }]])->{relationships} };
        $rel->{properties} = $updated_rel->{properties};
        return;
    }
    
    # incoming/outgoing/undirected hash refs are {min_depth, max_depth, type,
    # namespace, label, properties}, where the later 3 are result node specs and
    # with depths defaulting to 1 and others defaulting to undef; none supplied
    # defaults to undirected {min_depth => 1, max_depth => 1}. There is also a
    # left/rightmost boolean that combined with max_depth over 1 will return
    # just the left (for incoming) or right (for outgoing) most node in the path.
    # This returns hash of nodes and relationships for use by frontends;
    # related_nodes() calls this and returns just a list of nodes
    sub related {
        my ($self, $start_node, $undirected, $incoming, $outgoing, $result_nodes_only, $return_history_nodes) = @_;
        $self->throw("undirected is mutually exclusive of in/outgoing") if $undirected && ($outgoing || $incoming);
        if (!$outgoing && !$incoming && !$undirected) {
            $undirected = { min_depth => 1, max_depth => 1 };
        }
        
        my $start_id = $self->node_id($start_node);
        
        # if we have a max_depth > 6 and a node spec and no type, we're going to
        # repeat the query while incrementing max_depth by 1, starting from
        # $min_depth, until we get any nodes. This could be orders of magnitude
        # faster than doing 1..20. (We can't cope with both incoming and
        # outgoing have a max_depth > 6, because we can't easily know which side
        # any returned nodes were for)
        my $limit_depth     = 0;
        my $depth_increment = 0;
        while (1) {
            my ($cypher, $props, $low_depth, $high_depth);
            
            if ($undirected) {
                my ($result_node_spec, $properties, $type, $min_depth, $max_depth) = $self->_related_nodes_hashref_parse($undirected, 'param');
                my $return = $result_nodes_only ? 'DISTINCT u' : 'p';
                if ($max_depth > 6 && $result_node_spec && !$type) {
                    $low_depth  = $min_depth;
                    $high_depth = $max_depth;
                }
                if ($limit_depth) {
                    $min_depth = $limit_depth;
                    $max_depth = $limit_depth;
                }
                my $depth = ($min_depth == 1 && $max_depth == 1) ? '' : "*$min_depth..$max_depth";
                $cypher = "MATCH p = (start)-[$type$depth]-(u$result_node_spec) WHERE id(start) = $start_id RETURN $return";
                $props = { 'param' => $properties };
            }
            else {
                my (%all_properties, @return);
                my $most = '';
                my $left = '';
                if ($incoming) {
                    my ($result_node_spec, $properties, $type, $min_depth, $max_depth) = $self->_related_nodes_hashref_parse($incoming, 'left');
                    if ($max_depth > 6 && $result_node_spec && !$type) {
                        $low_depth  = $min_depth;
                        $high_depth = $max_depth;
                    }
                    if ($limit_depth) {
                        $min_depth = $limit_depth;
                        $max_depth = $limit_depth;
                    }
                    my $depth = ($min_depth == 1 && $max_depth == 1) ? '' : "*$min_depth..$max_depth";
                    $left = "(l$result_node_spec)-[$type$depth]->";
                    push(@return, 'l');
                    $all_properties{left} = $properties if $properties;
                    if ($incoming->{leftmost} && $max_depth > 1) {
                        $most = "AND NOT (l)<-[$type]-($result_node_spec) ";
                    }
                }
                my $right = '';
                if ($outgoing) {
                    my ($result_node_spec, $properties, $type, $min_depth, $max_depth) = $self->_related_nodes_hashref_parse($outgoing, 'right');
                    if ($max_depth > 6 && $result_node_spec && !$type) {
                        if ($high_depth) {
                            # can't cope with both incoming and outgoing having
                            # high max_depth, will behave as if neither did
                            undef $high_depth;
                        }
                        else {
                            $low_depth  = $min_depth;
                            $high_depth = $max_depth;
                        }
                    }
                    if ($limit_depth) {
                        $min_depth = $limit_depth;
                        $max_depth = $limit_depth;
                    }
                    my $depth = ($min_depth == 1 && $max_depth == 1) ? '' : "*$min_depth..$max_depth";
                    $right = "-[$type$depth]->(r$result_node_spec)";
                    push(@return, 'r');
                    $all_properties{right} = $properties if $properties;
                    if ($outgoing->{rightmost} && $max_depth > 1) {
                        $most = "AND NOT (r)-[$type]->($result_node_spec) ";
                    }
                }
                
                my $return;
                if ($result_nodes_only) {
                    $return = 'DISTINCT ' . join(', ', @return);
                }
                else {
                    $return = 'p';
                }
                
                $cypher = "MATCH p = $left(start)$right where id(start) = $start_id ${most}RETURN $return";
                $props = keys %all_properties ? \%all_properties : undef;
            }
            
            if ($high_depth) {
                if ($limit_depth) {
                    my $hash = $self->_run_cypher([[$cypher, $props ? ($props) : ()]], $return_history_nodes ? ({ return_history_nodes => 1 }) : ());
                    if ($hash && defined $hash->{nodes} && @{ $hash->{nodes} }) {
                        return $hash;
                    }
                }
                
                $limit_depth = $low_depth + $depth_increment++;
                last if $limit_depth > $high_depth;
                redo;
            }
            else {
                return $self->_run_cypher([[$cypher, $props ? ($props) : ()]], $return_history_nodes ? ({ return_history_nodes => 1 }) : ());
            }
        }
    }
    
    sub _related_nodes_hashref_parse {
        my ($self, $hashref, $param_key) = @_;
        my $type = $hashref->{type} ? ':' . join('|', map { "`$_`" } split(/\|/, $hashref->{type})) : '';
        my $min_depth = defined $hashref->{min_depth} ? $hashref->{min_depth} : 1;
        my $max_depth = $hashref->{max_depth} || $min_depth || 1;
        my $result_node_spec = '';
        if ($hashref->{namespace} && $hashref->{label}) {
            my ($labels, $param_map) = $self->_labels_and_param_map($hashref->{namespace}, $hashref->{label}, $hashref->{properties}, $param_key);
            $result_node_spec = ':' . $labels . $param_map;
        }
        return ($result_node_spec, $hashref->{properties}, $type, $min_depth, $max_depth);
    }
    
    method related_nodes (HashRef|Object $start_node!, Bool :$return_history_nodes = 0, HashRef :$outgoing?, HashRef :$incoming?, HashRef :$undirected?) {
        my $start_id = $start_node->{id};
        my $hash = $self->related($start_node, $undirected, $incoming, $outgoing, 1, $return_history_nodes);
        return unless ($hash && defined $hash->{nodes});
        return grep { $_->{id} != $start_id } @{ $hash->{nodes} };
    }
    
    # direction is 'incoming' or 'outgoing'. Undef means undirected
    # properties is [['key', 'value', 0], [ ... ]], where third value is a booleon which if true means the value is treated as a regex
    method closest_nodes_with_label (HashRef|Object $start_node!, Str $namespace!, Str $label!, Str :$direction?, ArrayRef[ArrayRef] :$properties?, Int :$depth = 100, Bool :$all = 0) {
        my $start_id = $start_node->{id};
        my $dir      = $direction ? "&direction=$direction" : '';
        my $prop     = '';
        if ($properties) {
            my @props;
            foreach my $prop_array (@$properties) {
                my ($key, $val, $regex) = @$prop_array;
                my $type = $regex ? 'regex' : 'literal';
                $val = uri_escape($val);
                push(@props, "$type\%40_\%40$key\%40_\%40$val");
            }
            my $props = join('%40%40%40', @props);
            $prop = "&properties=$props";
        }
        
        return $self->_call_vrpipe_neo4j_plugin_and_parse("/closest/$global_label\%7C$namespace\%7C$label/to/$start_id?depth=$depth&all=$all$dir$prop", namespace => $namespace, label => $label);
    }
    
    method _call_vrpipe_neo4j_plugin_and_parse (Str $path!, Str :$namespace!, Str :$label?) {
        my $data = $self->_call_vrpipe_neo4j_plugin($path);
        
        my @nodes;
        if ($data && ref($data) eq 'HASH') {
            while (my ($nid, $props) = each %{$data}) {
                my $this_label = delete $props->{neo4j_label};
                $this_label ||= $label;
                $this_label || $self->throw("No label supplied, and none returned by the plugin");
                push(@nodes, { id => int($nid), namespace => $namespace, label => $this_label, properties => $props });
            }
        }
        
        if (@nodes) {
            if (wantarray) {
                return @nodes;
            }
            return $nodes[0];
        }
        return;
    }
    
    sub _call_vrpipe_neo4j_plugin {
        # this plugin needs to be installed in Neo4J first:
        # https://github.com/VertebrateResequencing/vrpipe_neo4j_plugin
        my ($self, $path) = @_;
        my $pluginurl = "$url/v1/service$path";
        
        my $data;
        foreach my $try_num (1 .. 20) {
            eval {
                my $tx = $ua->get($pluginurl => $ua_headers);
                my $res = $tx->success;
                unless ($res) {
                    my $err = $tx->error;
                    $self->throw("Failed to connect to '$pluginurl': [$err->{code}] $err->{message}");
                }
                $data = $json->decode($res->body);
            };
            if ($@) {
                if ($try_num == 20) {
                    die $@;
                }
                else {
                    sleep($try_num * 2);
                }
            }
            else {
                last;
            }
        }
        
        return $data;
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
        my $str = shift;
        return unless $str;
        return $json->decode($str);
    }
    
    method date_to_epoch (Str $dstr) {
        # allow copy/pastes of stringified DateTimes
        $dstr =~ s/(\d)T(\d)/$1 $2/;
        
        my $parser = DateTime::Format::Natural->new;
        my $dt     = $parser->parse_datetime($dstr);
        if ($parser->success) {
            return $dt->epoch;
        }
        else {
            $self->throw($parser->error);
        }
    }
}

1;
