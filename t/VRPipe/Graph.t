#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Parallel::ForkManager;
use Path::Class;

BEGIN {
    use Test::Most tests => 33;
    use VRPipeTest;
    use_ok('VRPipe::Persistent::Graph');
}

ok my $graph = VRPipe::Persistent::Graph->new, 'able to create a Graph object';
isa_ok($graph, 'VRPipe::Persistent::Graph');

ok $graph->drop_database, 'we were able to drop the test database';

ok $graph->add_schema(namespace => 'VRTrack', label => 'Sample', unique => [qw(sanger_id)], indexed => [qw(uuid public_name)]), 'add_schema worked';
is $graph->add_schema(namespace => 'VRTrack', label => 'Sample', unique => [qw(sanger_id)], indexed => [qw(uuid public_name)]), 0, 'add_schema a second time with same args does nothing';
throws_ok { $graph->add_schema(namespace => 'OtherNS', label => 'Sam|ple', unique => [qw(name)]) } qr/neither namespace or label may contain the \| character/, 'add_schema throws if you try a label with the pipe character in it';
throws_ok { $graph->add_schema(namespace => 'OtherNS', label => 'Sample', unique => ['nam|e']) } qr/parameter may not contain the \| character/, 'add_schema throws if you try a parameter with the pipe character in it';
ok $graph->add_schema(namespace => 'OtherNS', label => 'Sample', unique => [qw(name)]), 'add_schema worked for the same label in a different namespace';
my ($unique, $indexed) = $graph->get_schema(namespace => 'VRTrack', label => 'Sample');
is_deeply [$unique, $indexed], [[qw(sanger_id)], [qw(uuid public_name)]], 'get_schema returned the unique and indexed fields from its internal cache';
$VRPipe::Persistent::Graph::schemas = {};
($unique, $indexed) = $graph->get_schema(namespace => 'VRTrack', label => 'Sample');
is_deeply [$unique, $indexed], [[qw(sanger_id)], [qw(uuid public_name)]], 'get_schema returned the unique and indexed fields from the database';

ok my $node = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' }), 'was able to add a node';
isa_ok($node, 'REST::Neo4p::Node');
my $orig_node_id = $node->id;
my $sanger1 = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' });
is $sanger1->id, $orig_node_id, 'calling add_node twice with the same args returns the same node';
throws_ok { $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', public_name => 'public2' }) } qr/already exists with label vdt.+\|VRTrack\|Sample/, 'add_node throws when used twice with the same unique arg but different other arg';
ok $graph->delete_node($sanger1), 'delete_node() worked';
$sanger1 = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' });
is $sanger1->id, $orig_node_id, 'calling add_node again with the same args after deleting returns a new node';
throws_ok { $graph->add_node(namespace => 'VRTrack', label => 'Sampl', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' }) } qr/You must first create a schema for namespace/, 'add_node throws when supplied with a label for which there is no schema';

# get nodes by their properties (first add some more for testing purposes)
my $sanger2 = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger2', public_name => 'public1' });
my $sanger3 = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger3', public_name => 'public2' });
my @nodes = $graph->get_nodes(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1' });
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger1'], 'get_nodes returned expected node when searching on a unique value';
@nodes = $graph->get_nodes(namespace => 'VRTrack', label => 'Sample', properties => { public_name => 'public1' });
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger1', 'sanger2'], 'get_nodes returned expected nodes when searching on a non-unique value';

# relationship tests (first add even more nodes to test with)
ok $graph->add_schema(namespace => 'VRTrack', label => 'Individual', unique => [qw(name)]), 'add_schema worked with different args and no indexed';
$graph->add_schema(namespace => 'VRTrack', label => 'Library', unique => [qw(name)]);
$graph->add_schema(namespace => 'VRTrack', label => 'Lane',    unique => [qw(name)]);
$graph->add_schema(namespace => 'VRTrack', label => 'Study',   unique => [qw(name)]);
my $john = $graph->add_node(namespace => 'VRTrack', label => 'Individual', properties => { name => 'John' });
my $jane = $graph->add_node(namespace => 'VRTrack', label => 'Individual', properties => { name => 'Jane' });
ok $graph->relate($jane, $sanger3, 'sampled'), 'relate() worked';
$graph->relate($john, $sanger2, 'sampled');
$graph->relate($john, $sanger1, 'sampled');
my $study = $graph->add_node(namespace => 'VRTrack', label => 'Study', properties => { name => 'Study of Disease_xyz' });
$graph->relate($study, $john, 'has_participant');
$graph->relate($study, $jane, 'has_participant');
my @samples = ($sanger1, $sanger2, $sanger3);

foreach my $i (1 .. 3) {
    my $library = $graph->add_node(namespace => 'VRTrack', label => 'Library', properties => { name => "Library$i" });
    $graph->relate($samples[$i - 1], $library, 'prepared');
    my $lane = $graph->add_node(namespace => 'VRTrack', label => 'Lane', properties => { name => "Lane$i", $i == 2 ? (foo => 'bar') : () });
    $graph->relate($library, $lane, 'sequenced');
}

@nodes = $graph->related_nodes($john, namespace => 'VRTrack', label => 'Sample');
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger1', 'sanger2'], 'related_nodes returned all Sample nodes related to John';
@nodes = $graph->related_nodes($john, namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger2' });
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger2'], 'related_nodes can be limited by params';
@nodes = $graph->related_nodes($john, namespace => 'VRTrack', label => 'Lane');
is_deeply [sort map { $_->get_property('name') } @nodes], [], 'related_nodes returned nothing looking for Lanes related to John because the default depth is only 1';
@nodes = $graph->related_nodes($john, namespace => 'VRTrack', label => 'Lane', max_depth => 3);
is_deeply [sort map { $_->get_property('name') } @nodes], [qw(Lane1 Lane2)], 'related_nodes returned the lanes related to John given a depth of 3';
@nodes = $graph->related_nodes($john, namespace => 'VRTrack', label => 'Lane', max_depth => 3, direction => '>');
is_deeply [sort map { $_->get_property('name') } @nodes], [qw(Lane1 Lane2)], 'related_nodes still works given an outgoing direction';
@nodes = $graph->related_nodes($john, namespace => 'VRTrack', label => 'Lane', max_depth => 3, direction => '<');
is_deeply [sort map { $_->get_property('name') } @nodes], [], 'related_nodes returns no lanes given an incoming direction';
@nodes = $graph->related_nodes($john, namespace => 'VRTrack', label => 'Study', direction => '<');
is_deeply [sort map { $_->get_property('name') } @nodes], ['Study of Disease_xyz'], 'related_nodes works with an incoming direction to find the study';
@nodes = $graph->related_nodes($john, namespace => 'VRTrack', label => 'Lane', max_depth => 3, direction => '>', properties => { foo => 'bar' });
is_deeply [sort map { $_->get_property('name') } @nodes], [qw(Lane2)], 'related_nodes still works when using max_depth, direction and properties all at the same time';

# add some more nodes to test the visualisation
$study = $graph->add_node(namespace => 'VRTrack', label => 'Study', properties => { name => 'Study of Disease_abc' });
$graph->relate($study, $jane, 'has_participant');
$graph->add_schema(namespace => 'OtherNS', label => 'Workplace', unique => [qw(name)], indexed => []);
my $workplace = $graph->add_node(namespace => 'OtherNS', label => 'Workplace', properties => { name => 'Sanger' });
$graph->add_schema(namespace => 'OtherNS', label => 'Individual', unique => [qw(name)], indexed => []);
for (1 .. 100) {
    $node = $graph->add_node(namespace => 'OtherNS', label => 'Individual', properties => { name => 'Person' . $_ });
    $graph->relate($workplace, $node, 'has_employee');
}

# add nodes to test display of images, and test the uuid creation helper method
# and required schema option
my ($lane1) = $graph->get_nodes(namespace => 'VRTrack', label => 'Lane', properties => { name => 'Lane1' });
$graph->add_schema(namespace => 'VRPipe', label => 'StepResult', unique => [qw(uuid)]);
$graph->add_schema(namespace => 'VRPipe', label => 'Image', unique => [qw(path)], required => ['type']);
my $result = $graph->add_node(namespace => 'VRPipe', label => 'StepResult', properties => { uuid => $graph->create_uuid });
like $result->get_property('uuid'), qr/\w{8}-\w{4}-\w{4}-\w{4}-\w{12}/, 'create_uuid() method worked to populate a property';
throws_ok { $graph->add_node(namespace => 'VRPipe', label => 'StepResult', properties => { foo => 'bar' }) } qr/Parameter 'uuid' must be supplied/, 'add_node throws when a unique parameter is not supplied';
my $image = $graph->add_node(namespace => 'VRPipe', label => 'Image', properties => { path => file(qw(t data qcgraph.png))->absolute->stringify, type => 'png' });
throws_ok { $graph->add_node(namespace => 'VRPipe', label => 'Image', properties => { path => '/foo' }) } qr/Parameter 'type' must be supplied/, 'add_node throws when a specified required parameter is not supplied';

exit;
