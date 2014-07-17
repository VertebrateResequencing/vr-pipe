
=head1 NAME

VRPipe::SchemaLabelRole - a role for auto-generated schema label classes

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This role is used by classes auto-generated during VRPipe::Schema->create().
See Schema docs for more information.

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

role VRPipe::SchemaLabelRole {
    use VRPipe::Persistent::Graph;
    my $graph = VRPipe::Persistent::Graph->new();
    
    method _get_setter (Str $property!, Str $new_value?) {
        if (defined $new_value) {
            $graph->node_add_properties($self, { $property => $new_value });
        }
        else {
            return $graph->node_property($self, $property);
        }
    }
    
    method node_id {
        return $self->{id};
    }
    
    method properties (Bool :$flatten_parents = 0) {
        return $graph->node_properties($self, flatten_parents => $flatten_parents);
    }
    
    method parent_property (Str $property) {
        my $properties = $self->properties(flatten_parents => 1);
        return $properties->{$property} if exists $properties->{$property};
        return;
    }
    
    method add_properties (HashRef $properties!) {
        $graph->node_add_properties($self, $properties);
    }
    
    method update_from_db {
        my $hash = $graph->get_node_by_id($self->node_id);
        $self->{properties} = $hash->{properties};
    }
    
    method related (HashRef :$outgoing?, HashRef :$incoming?, HashRef :$undirected?) {
        my @nodes = $graph->related_nodes($self, $outgoing ? (outgoing => $outgoing) : (), $incoming ? (incoming => $incoming) : (), $undirected ? (undirected => $undirected) : ());
        
        # bless them into the correct classes
        foreach my $node (@nodes) {
            bless $node, 'VRPipe::Schema::' . $node->{namespace} . '::' . $node->{label};
        }
        
        return @nodes;
    }
    
    method relate_to (HashRef|Object $node!, $type!) {
        $graph->relate($self, $node, type => $type);
    }
}

1;
