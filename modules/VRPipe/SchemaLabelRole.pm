
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
    use VRPipe::Persistent::InMemory;
    my $im    = VRPipe::Persistent::InMemory->new();
    my $graph = VRPipe::Persistent::Graph->new();
    my $pwh;
    
    has 'valid_properties' => (
        traits  => ['Hash'],
        is      => 'ro',
        isa     => 'HashRef[Str]',
        builder => '_build_valid_properties',
        lazy    => 1,
        handles => {
            property_is_valid => 'exists',
        },
    );
    
    method _build_valid_properties {
        return { map { $_ => 1 } ($self->unique_properties, $self->required_properties, $self->other_properties) };
    }
    
    method _get_setter (Str $property!, Str $new_value?) {
        # this should only be called by the auto-generated accessor methods, in
        # which case we don't need to check that $property is valid
        
        if (defined $new_value) {
            my %uniques = map { $_ => 1 } $self->unique_properties();
            if (exists $uniques{$property}) {
                my $class = $self->class();
                $self->throw("Property '$property' is unique for schema $class and can't be changed");
            }
            
            $self->block_until_locked();
            $graph->node_add_properties($self, { $property => $new_value });
            $self->_maintain_property_history(0);
            $self->unlock();
        }
        else {
            return $graph->node_property($self, $property);
        }
    }
    
    method node_id {
        return $self->{id};
    }
    
    method block_until_locked {
        my $lock_key = 'Schema.Node.' . $self->node_id();
        $im->block_until_locked($lock_key);
    }
    
    method unlock {
        my $lock_key = 'Schema.Node.' . $self->node_id();
        $im->unlock($lock_key);
    }
    
    method own_lock {
        my $lock_key = 'Schema.Node.' . $self->node_id();
        return $im->locked($lock_key, by_me => 1);
    }
    
    method properties (Bool :$flatten_parents = 0) {
        return $graph->node_properties($self, flatten_parents => $flatten_parents);
    }
    
    method parent_property (Str $property) {
        my $properties = $self->properties(flatten_parents => 1);
        return $properties->{$property} if exists $properties->{$property};
        return;
    }
    
    method add_properties (HashRef $properties!, Bool :$replace = 0) {
        my $allows_anything = $self->allows_anything();
        my $valid_props     = $self->valid_properties();
        my %uniques         = map { $_ => 1 } $self->unique_properties();
        my $class           = $self->class();
        foreach my $prop (keys %$properties) {
            unless ($allows_anything || exists $valid_props->{$prop}) {
                $self->throw("Property '$prop' supplied, but that isn't defined in the schema for $class");
            }
            if (exists $uniques{$prop}) {
                $self->throw("Property '$prop' supplied, but that's unique for schema $class and can't be changed");
            }
        }
        
        $self->block_until_locked();
        my $graph_method = $replace ? 'node_set_properties' : 'node_add_properties';
        
        if ($replace) {
            # add unique and missing required properties back in to what we'll
            # set
            my %props_to_set = map { $_ => $properties->{$_} } keys %$properties;
            foreach my $prop (keys %uniques, $self->required_properties()) {
                next if exists $props_to_set{$prop};
                $props_to_set{$prop} = $self->$prop();
            }
            $properties = \%props_to_set;
        }
        
        $graph->$graph_method($self, $properties);
        $self->_maintain_property_history($replace);
        $self->unlock();
    }
    
    method _maintain_property_history (Bool $replace) {
        return unless $self->_keep_history;
        
        my %uniques = map { $_ => 1 } $self->unique_properties();
        my $history_props;
        while (my ($key, $val) = each %{ $self->{properties} }) {
            unless (exists $uniques{$key}) {
                $history_props->{$key} = $val;
            }
        }
        return unless $history_props;
        
        $pwh ||= VRPipe::Schema->create('PropertiesWithHistory');
        my $changed_properties = $pwh->_add_or_update_properties($graph, $self, $history_props, replace => $replace);
        $self->{changed_properties} = $changed_properties;
    }
    
    method changed {
        if (defined $self->{changed_properties}) {
            return $self->{changed_properties};
        }
        return;
    }
    
    # returns a list of hashrefs with keys group_uuid, timestamp, properties,
    # where properties is another hashref of the actual properties at that
    # timestamp (excluding unique properties, which don't change)
    method property_history (Str $property?) {
        $pwh ||= VRPipe::Schema->create('PropertiesWithHistory');
        return $pwh->_get_property_history($graph, $self, $property ? ($property) : ());
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
    
    method relate_to (HashRef|Object $node!, $type!, Bool :$selfish = 0, Bool :$replace = 0) {
        $graph->relate($self, $node, type => $type, selfish => $selfish, replace => $replace);
    }
}

1;
