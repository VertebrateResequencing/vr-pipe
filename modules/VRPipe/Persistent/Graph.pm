
=head1 NAME

VRPipe::Persistent::Graph - interface to a graph database

=head1 SYNOPSIS
    
[...]

=head1 DESCRIPTION

For schema-less store of connected data we use a graph database; Neo4J in this
case.

This is essentially a wrapper around REST::Neo4p, providing functions that can
be used to store and retrieve information about things.

Things (must) have labels and properties, and can have dynamically-applied
schema restrictions to provide uniqueness and indexes.

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
    use REST::Neo4p;
    
    our $vrp_config = VRPipe::Config->new();
    our ($neo4p, $global_label);
    
    sub BUILD {
        my $self = shift;
        
        unless ($neo4p) {
            my $url = $vrp_config->neo4j_server_url();
            $neo4p = REST::Neo4p->connect($url);
            
            my $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
            if ($deployment eq 'production') {
                $global_label = "deployment|production";
            }
            else {
                my $user = getlogin || getpwuid($<);
                $global_label = "deployment|testing|$user";
            }
        }
    }
    
    sub _run_cypher {
        my ($self, $cypher, $params) = @_;
        my $query = REST::Neo4p::Query->new($cypher, $params ? ({ param => $params }) : ());
        $query->execute;
        if ($query->err) {
            $self->throw("neo4j cypher query failed with code " . $query->err . ": " . $query->errstr);
        }
        else {
            my @results;
            while (my $row = $query->fetch) {
                push(@results, @$row > 1 ? $row : $row->[0]);
                #*** when we ask for a node neo4j will return the full details of the node, but it seems like neo4p only creates an object consisting of the id? maybe we should avoid object construction entirely and only return ID(node)?
            }
            return @results;
        }
    }
    
    method drop_database {
        $self->throw("drop_database() can only be used when testing") unless $global_label =~ /^deployment\|testing/;
        $self->_run_cypher("MATCH (n:`$global_label`) OPTIONAL MATCH (n:`$global_label`)-[r]-() DELETE n,r");
        return 1;
    }
    
    method add_schema (Str :$label!, ArrayRef[Str] :$unique!, ArrayRef[Str] :$indexed?) {
        # have we already done this?
        my $schema_labels = qq[`$global_label`:`VRPipeInternals`:`Schema|$label`];
        my ($done) = $self->_run_cypher("MATCH (n:$schema_labels) RETURN n");
        unless ($done) {
            # set constraints (which also adds an index on the constraint)
            foreach my $field (@$unique) {
                $self->_run_cypher("CREATE CONSTRAINT ON (n:`$label`) ASSERT n.$field IS UNIQUE");
            }
            
            # add indexes
            foreach my $field (@{ $indexed || [] }) {
                $self->_run_cypher("CREATE INDEX ON :`$label`($field)");
            }
            
            # record that we've done this
            $self->_run_cypher("CREATE (:$schema_labels)");
            
            return 1;
        }
        return 0;
    }
    
    sub _labels_and_param_map {
        my ($self, $label, $params) = @_;
        my $labels = "`$global_label`:`VRPipe`:`$label`";
        my $param_map = $params ? ' { ' . join(', ', map { "$_: {param}.$_" } sort keys %$params) . ' }' : '';
        return ($labels, $param_map);
    }
    
    method add_node (Str :$label!, HashRef :$properties!) {
        my ($labels, $param_map) = $self->_labels_and_param_map($label, $properties);
        
        if (defined wantarray()) {
            my ($node) = $self->_run_cypher("MERGE (n:$labels$param_map) RETURN n", $properties);
            return $node;
        }
        else {
            $self->_run_cypher("MERGE (:$labels$param_map)", $properties);
        }
    
    }
    
    method get_nodes (Str :$label!, HashRef :$properties!) {
        my ($labels, $param_map) = $self->_labels_and_param_map($label, $properties);
        return $self->_run_cypher("MATCH (n:$labels$param_map) RETURN n", $properties);
    }
    
    method relate (Object $start_node!, Object $end_node!, Str $relationship!) {
        #*** is there a benefit to constructing our own cypher query here instead?
        return $start_node->relate_to($end_node, $relationship);
    }
    
    method related_nodes (Object $start_node!, Str :$label!, HashRef :$properties?, Str :$relationship?, Str :$direction?, Int :$min_depth = 1, Int :$max_depth = 1) {
        my ($labels, $param_map) = $self->_labels_and_param_map($label, $properties);
        my $type = $relationship ? ":`$relationship`" : '';
        $direction ||= '';
        my $leftward  = $direction eq '<' ? '<' : '';
        my $rightward = $direction eq '>' ? '>' : '';
        return $self->_run_cypher("START start=node(" . $start_node->id . ") MATCH (start)$leftward-[$type*$min_depth..$max_depth]-$rightward(n:$labels$param_map) RETURN n", $properties);
    }
}

1;
