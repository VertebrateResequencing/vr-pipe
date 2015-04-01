
=head1 NAME

VRPipe::Schema::PropertiesWithHistory - schemas for storing node properties

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This simple Schema lets VRPipe store node properties while keeping their
history. Internal-use only by VRPipe::SchemaRole.

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

class VRPipe::Schema::PropertiesWithHistory with VRPipe::SchemaRole {
    use Digest::MD5 qw(md5_hex);
    use VRPipe::Persistent::InMemory;
    
    method schema_definitions {
        return [{
                label    => 'PropertyGroup',
                unique   => [qw(uuid)],
                required => [qw(timestamp)]
            },
            {
                label   => 'Property',
                unique  => [qw(keyval_md5)],
                indexed => [qw(key value)]
            }
        ];
    }
    
    method _nodes_to_properties (ArrayRef $nodes) {
        my $property_group;
        my $properties = {};
        my %property_node_details;
        foreach my $n (@$nodes) {
            if ($n->{namespace} eq 'PropertiesWithHistory') {
                if ($n->{label} eq 'Property') {
                    my $p   = $n->{properties};
                    my $key = $p->{key};
                    my $val = $p->{value};
                    if (exists $properties->{$key}) {
                        my $current = $properties->{$key};
                        if (ref($current)) {
                            push(@$current, $val);
                        }
                        else {
                            $properties->{$key} = [$current, $val];
                        }
                    }
                    else {
                        $properties->{$key} = $val;
                    }
                    
                    my $md5 = $p->{keyval_md5};
                    $property_node_details{$key}->{$md5} = { keyval_md5 => $md5, key => $key, value => $val };
                }
                else {
                    $property_group = $n;
                }
            }
        }
        
        return ($property_group, \%property_node_details, $properties);
    }
    
    method _add_or_update_properties ($graph, $node, HashRef $new_props, Bool :$replace = 0) {
        my $locked_by_me = $node->own_lock();
        $node->block_until_locked unless $locked_by_me;
        
        # get the current property group and nodes (if any)
        my $node_id = $node->node_id();
        my $data = $graph->_run_cypher([["MATCH (n)-[:current_properties*1]->(c)-[:property*1]->(p) WHERE id(n) = $node_id RETURN c,p"]], { return_history_nodes => 1 });
        my ($current_property_group, $current_property_node_details, $current_props) = $self->_nodes_to_properties($data->{nodes});
        
        # calculate the property node details we'll need for the new properties
        my %needed_prop_node_details;
        while (my ($key, $val) = each %$new_props) {
            my @new_vals = ref($val) ? (sort @$val) : ($val);
            foreach my $v (@new_vals) {
                my $md5 = md5_hex($key . '::' . $v);
                $needed_prop_node_details{$key}->{$md5} = { keyval_md5 => $md5, key => $key, value => $v };
            }
        }
        
        # note the critical differences that would require us to form a new
        # property group
        my ($changed, @final_group);
        while (my ($key, $hash) = each %needed_prop_node_details) {
            push(@final_group, values %$hash);
            
            if ($current_property_group && exists $current_property_node_details->{$key}) {
                my @current = sort keys %{ $current_property_node_details->{$key} };
                my @new     = sort keys $needed_prop_node_details{$key};
                if ("@current" ne "@new") {
                    $changed->{$key} = [$current_props->{$key}, $new_props->{$key}];
                }
            }
            else {
                $changed->{$key} = [undef, $new_props->{$key}];
            }
        }
        
        # include in the final group current properties not mentioned in
        # new, or just note that we don't have them any more if in replace mode
        if ($current_property_group) {
            while (my ($key, $hash) = each %$current_property_node_details) {
                next if exists $needed_prop_node_details{$key};
                if ($replace) {
                    $changed->{$key} = [$current_props->{$key}, undef];
                }
                else {
                    push(@final_group, values %$hash);
                }
            }
        }
        
        unless ($changed) {
            $node->unlock unless $locked_by_me;
            return;
        }
        
        # make a new property group and attach to (potentially new) property
        # nodes in the final grouping
        my $uuid               = $graph->create_uuid;
        my $time               = time();
        my $new_property_group = $self->add('PropertyGroup', { uuid => $uuid, timestamp => $time }, $current_property_group ? (outgoing => { node => $current_property_group, type => 'previous_properties' }) : ());
        $node->relate_to($new_property_group, 'current_properties', replace => 1);
        
        # lock before possibly creating properties to avoid another process
        # deleting them because it is deleting some other node that shares a
        # some properties; this should be a very rare event so its ok to lock on
        # a very general key
        my $im       = VRPipe::Persistent::InMemory->new();
        my $lock_key = 'graph.propertieswithhistory.updating';
        $im->block_until_locked($lock_key);
        $self->add('Property', \@final_group, incoming => { node => $new_property_group, type => 'property' });
        $im->unlock($lock_key);
        
        $node->unlock unless $locked_by_me;
        return $changed;
    }
    
    method _get_property_history ($graph, $node, Str $property?) {
        my $node_id = $node->node_id();
        #*** arbitrarily limited depth of history to 500; we should do paging or
        #    make it an option or something
        my $property_where = $property ? " AND p.key = '$property'" : '';
        my $data = $graph->_run_cypher([["MATCH (n)-[:current_properties*1]->(c)-[r:property*1]->(p) WHERE id(n) = $node_id$property_where OPTIONAL MATCH (c)-[:previous_properties*1..500]->()-[s:property*1]->(q) RETURN c,r,p,s,q"]], { return_history_nodes => 1 });
        
        # based on the relationships we can group the properties with their
        # propertygroup (properties may belong to more than 1 group)
        my %nodes;
        foreach my $node (@{ $data->{nodes} }) {
            $nodes{ $node->{id} } = $node;
        }
        
        my %groups;
        foreach my $rel (@{ $data->{relationships} }) {
            push(@{ $groups{ $rel->{startNode} } }, $nodes{ $rel->{endNode} });
        }
        
        # sort and generate desired output
        my @history;
        foreach my $group_node_id (sort { $b <=> $a } keys %groups) {
            my $property_nodes            = $groups{$group_node_id};
            my $property_group_properties = $nodes{$group_node_id}->{properties};
            my (undef, undef, $properties) = $self->_nodes_to_properties($property_nodes);
            push(@history, { properties => $properties, timestamp => $property_group_properties->{timestamp}, group_uuid => $property_group_properties->{uuid} });
        }
        
        return @history;
    }
}

1;
