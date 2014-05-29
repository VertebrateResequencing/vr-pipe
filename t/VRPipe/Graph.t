#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Parallel::ForkManager;

BEGIN {
    use Test::Most tests => 22;
    use VRPipeTest;
    use_ok('VRPipe::Persistent::Graph');
}

ok my $graph = VRPipe::Persistent::Graph->new, 'able to create a Graph object';
isa_ok($graph, 'VRPipe::Persistent::Graph');

ok $graph->drop_database, 'we were able to drop the test database';

ok $graph->add_schema(label => 'Sample', unique => [qw(sanger_id)], indexed => [qw(public_name)]), 'add_schema worked';
is $graph->add_schema(label => 'Sample', unique => [qw(sanger_id)], indexed => [qw(public_name)]), 0, 'add_schema a second time with same args does nothing';

ok my $node = $graph->add_node(label => 'Sample', properties => { sanger_id => 'sanger1', public_name => 'public1' }), 'was able to add a node';
isa_ok($node, 'REST::Neo4p::Node');
my $orig_node_id = $node->id;
my $sanger1 = $graph->add_node(label => 'Sample', properties => { sanger_id => 'sanger1', public_name => 'public1' });
is $sanger1->id, $orig_node_id, 'calling add_node twice with the same args returns the same node';
throws_ok { $graph->add_node(label => 'Sample', properties => { sanger_id => 'sanger1', public_name => 'public2' }) } qr/already exists with label Sample/, 'add_node throws when used twice with the same unique arg but different other arg';

# get nodes by their properties (first add some more for testing purposes)
my $sanger2 = $graph->add_node(label => 'Sample', properties => { sanger_id => 'sanger2', public_name => 'public1' });
my $sanger3 = $graph->add_node(label => 'Sample', properties => { sanger_id => 'sanger3', public_name => 'public2' });
my @nodes = $graph->get_nodes(label => 'Sample', properties => { sanger_id => 'sanger1' });
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger1'], 'get_nodes returned expected node when searching on a unique value';
@nodes = $graph->get_nodes(label => 'Sample', properties => { public_name => 'public1' });
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger1', 'sanger2'], 'get_nodes returned expected nodes when searching on a non-unique value';

# relationship tests (first add even more nodes to test with)
ok $graph->add_schema(label => 'Individual', unique => [qw(name)]), 'add_schema worked with different args and no indexed';
$graph->add_schema(label => 'Library', unique => [qw(name)]);
$graph->add_schema(label => 'Lane',    unique => [qw(name)]);
$graph->add_schema(label => 'Study',   unique => [qw(name)]);
my $john = $graph->add_node(label => 'Individual', properties => { name => 'John' });
my $jane = $graph->add_node(label => 'Individual', properties => { name => 'Jane' });
ok $graph->relate($sanger3, $jane, 'extracted_from'), 'relate() worked';
$graph->relate($sanger2, $john, 'extracted_from');
$graph->relate($sanger1, $john, 'extracted_from');
my $study = $graph->add_node(label => 'Study', properties => { name => 'Study of Disease_xyz' });
$graph->relate($john, $study, 'member_of');
$graph->relate($jane, $study, 'member_of');
my @samples = ($sanger1, $sanger2, $sanger3);

foreach my $i (1 .. 3) {
    my $library = $graph->add_node(label => 'Library', properties => { name => "Library$i" });
    $graph->relate($library, $samples[$i - 1], 'made_from');
    my $lane = $graph->add_node(label => 'Lane', properties => { name => "Lane$i", $i == 2 ? (foo => 'bar') : () });
    $graph->relate($lane, $library, 'sequenced_from');
}

@nodes = $graph->related_nodes($john, label => 'Sample');
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger1', 'sanger2'], 'related_nodes returned all Sample nodes related to John';
@nodes = $graph->related_nodes($john, label => 'Sample', properties => { sanger_id => 'sanger2' });
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger2'], 'related_nodes can be limited by params';
@nodes = $graph->related_nodes($john, label => 'Lane');
is_deeply [sort map { $_->get_property('name') } @nodes], [], 'related_nodes returned nothing looking for Lanes related to John because the default depth is only 1';
@nodes = $graph->related_nodes($john, label => 'Lane', max_depth => 3);
is_deeply [sort map { $_->get_property('name') } @nodes], [qw(Lane1 Lane2)], 'related_nodes returned the lanes related to John given a depth of 3';
@nodes = $graph->related_nodes($john, label => 'Lane', max_depth => 3, direction => '<');
is_deeply [sort map { $_->get_property('name') } @nodes], [qw(Lane1 Lane2)], 'related_nodes still works given an incoming direction';
@nodes = $graph->related_nodes($john, label => 'Lane', max_depth => 3, direction => '>');
is_deeply [sort map { $_->get_property('name') } @nodes], [], 'related_nodes returns no lanes given an outgoing direction';
@nodes = $graph->related_nodes($john, label => 'Study', direction => '>');
is_deeply [sort map { $_->get_property('name') } @nodes], ['Study of Disease_xyz'], 'related_nodes works with an outgoing direction to find the study';
@nodes = $graph->related_nodes($john, label => 'Lane', max_depth => 3, direction => '<', properties => { foo => 'bar' });
is_deeply [sort map { $_->get_property('name') } @nodes], [qw(Lane2)], 'related_nodes still works when using max_depth, direction and properties all at the same time';

exit;

# # empty out the database
# my $cypher = "MATCH (n)
# OPTIONAL MATCH (n)-[r]-()
# DELETE n,r";
# execute($cypher);

# # try some stuff
# my $node_idx = REST::Neo4p::Index->new('node', 'named_index'); # old-style indexes, do not do this
# my $n1 = REST::Neo4p::Node->new( {name => 'Ferb'} );
# $node_idx->add_entry($n1, person => 'Ferb');
# my $n2 = REST::Neo4p::Node->new( {name => 'Phineas'} );
# $node_idx->add_entry($n2, person => 'Phineas');
# my $n3 = REST::Neo4p::Node->new( {name => 'Perry'} );
# $node_idx->add_entry($n3, animal => 'Perry');
# $n1->relate_to($n2, 'brother');
# $n3->relate_to($n1, 'pet');
# $n3->set_property({ species => 'Ornithorhynchus anatinus' });

# my ($result_node) = $node_idx->find_entries(person => 'Phineas');
# if ($result_node) {
#     $result_node->set_property({ species => 'Homo sapiens' });
# }
# $n2->set_property({ size => 'large' });

# # match
# $cypher = qq[MATCH (person { name: "Phineas" }) RETURN person;];
# ($result_node) = execute($cypher);
# $result_node->set_property({ colour => 'blue' });

# # create & merge
# $cypher = q[CREATE (n:Actor { name:"Tom Hanks" }) return n;];
# ($result_node) = execute($cypher);
# $result_node->set_property({ colour => 'red' });

# $cypher = q[CREATE (n:Actor { name:"Tom Hanks" }) return n;];
# ($result_node) = execute($cypher);
# $result_node->set_property({ dance => 'rumba' }) if $result_node;

# $cypher = q[CREATE CONSTRAINT ON (movie:Movie) ASSERT movie.title IS UNIQUE]; # this also creates an index
# # CREATE INDEX ON :Actor(name)
# execute($cypher);

# $cypher = q[MERGE (movie:Movie { title:"Fooballs", colour:"pink" }) ON CREATE SET movie.created = timestamp() ON MATCH set movie.matched = timestamp() return movie;];
# ($result_node) = execute($cypher); # creates
# $result_node->set_property({ got_fetch => '1' });
# ($result_node) = execute($cypher); # matches
# $result_node->set_property({ got_fetch => '2' });

# $cypher = q[MERGE (movie:Movie { title:"Fooballs", date:"2014" }) ON CREATE SET movie.created = timestamp() ON MATCH set movie.matched_again = timestamp() return movie;];
# ($result_node) = execute($cypher); # fails and does nothing due to uniqueness constraint on Movie.title
# $result_node->set_property({ got_fetch => '3' }) if $result_node;

# $cypher = q[CREATE (n1:Person { props1 }),(n2:Person { props2 }) RETURN n1,n2];
# (undef, $result_node) = execute($cypher, {props1 => { name => 'Andre', jacket => 'longsleeve' }, props2 => { name => 'Largo', jacket => 'shortsleeve' }});
# $result_node->set_property({ got_fetch => 'hoorah' }) if $result_node;

# $cypher = q[CREATE (n:Person { props }) RETURN n];
# (undef, $result_node) = execute($cypher, {props => [{ name => 'Joker', jacket => 'sleeveless' }, { name => 'Batman', jacket => 'kevlar' }]}); # only getting 1 result, don't know why (though both nodes get created)
# $result_node->set_property({ got_fetch => 'booyah' }) if $result_node; # doesn't happen :(

# $cypher = q[MERGE (n:Person { name: {param}.name, role: {param}.role }) RETURN n];
# ($result_node) = execute($cypher, { param => { name => 'Jim', role => 'ex-boss' }});
# $result_node->set_property({ got_merge_param => 'super' });

# $cypher = q[MATCH (a:Person { name: "Andre" }) RETURN a.jacket AS SomethingTotallyDifferent];
# my @results = execute($cypher); # [ [ 'longsleeve' ] ]
# #print Dumper(\@results), "\n";

# # set
# #*** should test if it's faster to use the set_property method, or to use cypher
# # set manually with params:
# $cypher = q[MATCH (n:Person { name: 'Andre' }) SET n.surname = { surname } RETURN n];
# execute($cypher, { surname => 'Giant' });

# $cypher = q[MATCH (n:Person { name: 'Joker' }) SET n = { props } RETURN n]; # removes all existing properties, sets the new ones
# execute($cypher, { props => { name => 'Joker', wiggle => 'cat' } });

# # delete
# $cypher = q[MATCH (n { name: 'Peter' }) DELETE n]; # delete a node
# execute($cypher); # no error if Peter didn't exist

# $cypher = "MATCH (n { name: 'Andre' })-[r]-() DELETE n, r";
# execute($cypher); # doesn't delete if Andre had no relationships

# #*** which is faster, getting using the method for a, or b?
# $cypher = "MATCH (a:Person { name:'Andre' }),(b:Person) WHERE b.name = 'Largo' CREATE (a)-[:KNOWS]->(b)";
# execute($cypher);

# $cypher = "MATCH (n { name: 'Andre' }) DELETE n";
# execute($cypher); # doesn't delete if Andre has a relationship

# $cypher = "MATCH (n { name: 'Andre' }) OPTIONAL MATCH (n)-[r]-() DELETE n, r";
# #execute($cypher); # sure-fire way to delete

# exit;

# sub execute {
#     my ($cypher, $params) = @_;
#     my $query = REST::Neo4p::Query->new($cypher, $params);
#     $query->execute;
#     if ($query->err) {
#         printf "status code: %d\n", $query->err;
#         printf "error message: %s\n", $query->errstr;
#         return;
#     }
#     else {
#         return @{$query->fetch || []};
#     }
# }
