
=head1 NAME

VRPipe::SchemaRole - a role that must be used by all Schemas

=head1 SYNOPSIS
    
    use VRPipe::Base;
    class VRPipe::Schema::MySchema with VRPipe::SchemaRole {
        method schema_definitions {
            return [
                {
                    label => 'Building',
                    unique => [qw(coordinates)],
                    indexed => [qw(name)],
                    required => [qw(name)],
                    optional => [qw(size)],
                }
            ];
        }
    }
    1;
    
    # then users can (see VRPipe::Schema docs for more):
    my $myschema = VRPipe::Schema->create('MySchema');
    my $building = $myschema->get('Building', { coordinates => 'xyz' });
    my $name = $building->name();
    $building->name('foo');

=head1 DESCRIPTION

Have a schema_definitions method that returns an array ref of hash refs, where
the hashes are the args you would supply to
VRPipe::Persistent::Graph->add_schema(), except for namespace, which will be
the class name. You can also add an additional 'optional' key to specify other
allowed attributes. Any attribute not specified somewhere in unique/indexed/
required/optional will not be get/settable via the auto-created
VRPipe::Schemas::MySchema::[label] class to which returned nodes will belong.

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

role VRPipe::SchemaRole {
    use VRPipe::Persistent::Graph;
    my $graph = VRPipe::Persistent::Graph->new();
    
    has 'schemas' => (
        is      => 'ro',
        isa     => 'ArrayRef',
        lazy    => 1,
        builder => 'schema_definitions'
    );
    
    has 'namespace' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_build_namespace'
    );
    
    has 'labels' => (
        traits  => ['Hash'],
        is      => 'ro',
        isa     => 'HashRef[ArrayRef[Str]]',
        default => sub { {} },
        handles => {
            _add_label       => 'set',
            valid_label      => 'exists',
            label_properties => 'get',
        },
    );
    
    method _build_namespace {
        my ($namespace) = ref($self) =~ /::([^:]+)$/;
        return $namespace;
    }
    
    method add_schemas {
        my $graph     = VRPipe::Persistent::Graph->new();
        my $namespace = $self->namespace;
        
        my $get_setter = sub {
        };
        
        foreach my $def (@{ $self->schemas }) {
            # create constraints and indexes in the database
            my $optional = delete $def->{optional};
            $graph->add_schema(%$def, namespace => $namespace);
            
            # store on ourselves what's valid according to this definition
            my $label = $def->{label};
            my @label_properties = (@{ $def->{unique} || [] }, @{ $def->{indexed} || [] }, @{ $def->{required} || [] }, @{ $optional || [] });
            $self->_add_label($label => \@label_properties);
            
            # create a class for this label
            my $methods = {};
            foreach my $property (@label_properties) {
                my $sub = sub {
                    shift->_get_setter($property, @_);
                };
                $methods->{ lc($property) } = $sub;
            }
            
            my $class = Moose::Meta::Class->create(
                'VRPipe::Schema::' . $namespace . '::' . $label,
                roles   => ['VRPipe::SchemaLabelRole'],
                methods => $methods
            );
        }
    }
    
    method graph {
        return $graph;
    }
    
    method _get_and_bless_nodes (Str $label!, Str $graph_method!, HashRef[Str]|ArrayRef[HashRef[Str]] $properties?, HashRef $extra_graph_args?) {
        my $namespace = $self->namespace;
        unless ($self->valid_label($label)) {
            $self->throw("'$label' isn't a valid label for schema $namespace");
        }
        
        my $props;
        if ($properties) {
            # check the supplied properties are allowed ($graph checks that required
            # ones are supplied)
            my %valid_props = map { $_ => 1 } @{ $self->label_properties($label) };
            $props = ref($properties) eq 'HASH' ? [$properties] : $properties;
            foreach my $prop_hash (@$props) {
                foreach my $prop (keys %$prop_hash) {
                    unless (exists $valid_props{$prop}) {
                        $self->throw("Property '$prop' supplied, but that isn't defined in the schema for ${namespace}::$label");
                    }
                }
            }
            
            if ($graph_method eq 'get_nodes') {
                $props = $properties;
            }
        }
        
        my @nodes = $graph->$graph_method(namespace => $namespace, label => $label, $props ? (properties => $props, ($graph_method eq 'add_nodes' ? (update => 1) : ())) : (), %{ $extra_graph_args || {} });
        
        # bless the nodes into the appropriate class
        foreach my $node (@nodes) {
            bless $node, 'VRPipe::Schema::' . $namespace . '::' . $label;
        }
        
        if (wantarray()) {
            return @nodes;
        }
        else {
            return $nodes[0];
        }
    }
    
    method add (Str $label!, HashRef[Str]|ArrayRef[HashRef[Str]] $properties!, HashRef :$incoming?, HashRef :$outgoing?) {
        return $self->_get_and_bless_nodes($label, 'add_nodes', $properties, { $incoming ? (incoming => $incoming) : (), $outgoing ? (outgoing => $outgoing) : () });
    }
    
    method get (Str $label!, HashRef $properties?) {
        return $self->_get_and_bless_nodes($label, 'get_nodes', $properties ? ($properties) : ());
    }
    
    method delete ($node) {
        $graph->delete_node($node);
    }
}

1;
