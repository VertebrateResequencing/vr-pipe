
=head1 NAME

VRPipe::Schema - apply schemas to the database and get/create/modify nodes

=head1 SYNOPSIS
    
    use VRPipe::Schema;
    
    my $schema = VRPipe::Schema->create($type);
    my $sample = $schema->add('Sample', { name => 'zyx' });
    my $library = $schema->add(
        'Library',
        { name => 'lkj' },
        incoming => { type => 'prepared', node => $sample }
    );
    $library = $schema->get('Library', { name => 'lkj' });
    
    my $name = $library->name();
    my $sample_name = $library->parent_properties('sample_name');
    $library->name('hgf');
    ($sample) = $library->related();
    
    $schema->delete($library);
    
    my $graph = $schema->graph;
    
    # to actually add a new schema to the database (or update a changed one):
    $schema = VRPipe::Schema->create($type, update_schemas_in_db => 1);

=head1 DESCRIPTION

A Schemas class holds a collection of related schemas. Calling create() for a
particular class will ensure that all those schemas are set in the database and
that you can now use VRPipe::Persistent::Graph methods to create and get nodes
under the namespace(s) and labels defined by the schemas.

create() returns an object that gives you convenience methods for adding,
getting and deleting nodes under one of the defined schemas. Unlike normal
graph nodes which are just plain hash references, the nodes returned by add()
and get() have methods on them named after each property of the label the node
belongs to.

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

package VRPipe::Schema;
use VRPipe::Base::AbstractFactory;

implementation_does qw/VRPipe::SchemaRole/;

sub create {
    my ($self, $type, %args) = @_;
    my $obj = $self->SUPER::create($type, {});
    $obj->add_schemas(%args);
    return $obj;
}

1;
