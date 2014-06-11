#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Parallel::ForkManager;

BEGIN {
    use Test::Most tests => 28;
    use VRPipeTest;
    use_ok('VRPipe::Persistent::Graph');
}

ok my $graph = VRPipe::Persistent::Graph->new, 'able to create a Graph object';
isa_ok($graph, 'VRPipe::Persistent::Graph');

ok $graph->drop_database, 'we were able to drop the test database';

ok $graph->add_schema(namespace => 'QCGrind', label => 'Sample', unique => [qw(sanger_id uuid)], indexed => [qw(public_name)]), 'add_schema worked';
is $graph->add_schema(namespace => 'QCGrind', label => 'Sample', unique => [qw(sanger_id uuid)], indexed => [qw(public_name)]), 0, 'add_schema a second time with same args does nothing';
throws_ok { $graph->add_schema(namespace => 'VRTrack', label => 'Sam|ple', unique => [qw(name)]) } qr/neither namespace or label may contain the \| character/, 'add_schema throws if you try a label with the pipe character in it';
throws_ok { $graph->add_schema(namespace => 'VRTrack', label => 'Sample', unique => ['nam|e']) } qr/parameter may not contain the \| character/, 'add_schema throws if you try a parameter with the pipe character in it';
ok $graph->add_schema(namespace => 'VRTrack', label => 'Sample', unique => [qw(name)]), 'add_schema worked for the same label in a different namespace';
my ($unique, $indexed) = $graph->get_schema(namespace => 'QCGrind', label => 'Sample');
is_deeply [$unique, $indexed], [[qw(sanger_id uuid)], ['public_name']], 'get_schema returned the unique and indexed fields from its internal cache';
$VRPipe::Persistent::Graph::schemas = {};
($unique, $indexed) = $graph->get_schema(namespace => 'QCGrind', label => 'Sample');
is_deeply [$unique, $indexed], [[qw(sanger_id uuid)], ['public_name']], 'get_schema returned the unique and indexed fields from the database';

ok my $node = $graph->add_node(namespace => 'QCGrind', label => 'Sample', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' }), 'was able to add a node';
isa_ok($node, 'REST::Neo4p::Node');
my $orig_node_id = $node->id;
my $sanger1 = $graph->add_node(namespace => 'QCGrind', label => 'Sample', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' });
is $sanger1->id, $orig_node_id, 'calling add_node twice with the same args returns the same node';
throws_ok { $graph->add_node(namespace => 'QCGrind', label => 'Sample', properties => { sanger_id => 'sanger1', public_name => 'public2' }) } qr/already exists with label QCGrind\|Sample/, 'add_node throws when used twice with the same unique arg but different other arg';
throws_ok { $graph->add_node(namespace => 'QCGrind', label => 'Sampl', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' }) } qr/You must first create a schema for namespace/, 'add_node throws when supplied with a label for which there is no schema';

# get nodes by their properties (first add some more for testing purposes)
my $sanger2 = $graph->add_node(namespace => 'QCGrind', label => 'Sample', properties => { sanger_id => 'sanger2', public_name => 'public1' });
my $sanger3 = $graph->add_node(namespace => 'QCGrind', label => 'Sample', properties => { sanger_id => 'sanger3', public_name => 'public2' });
my @nodes = $graph->get_nodes(namespace => 'QCGrind', label => 'Sample', properties => { sanger_id => 'sanger1' });
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger1'], 'get_nodes returned expected node when searching on a unique value';
@nodes = $graph->get_nodes(namespace => 'QCGrind', label => 'Sample', properties => { public_name => 'public1' });
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger1', 'sanger2'], 'get_nodes returned expected nodes when searching on a non-unique value';

# relationship tests (first add even more nodes to test with)
ok $graph->add_schema(namespace => 'QCGrind', label => 'Individual', unique => [qw(name)]), 'add_schema worked with different args and no indexed';
$graph->add_schema(namespace => 'QCGrind', label => 'Library', unique => [qw(name)]);
$graph->add_schema(namespace => 'QCGrind', label => 'Lane',    unique => [qw(name)]);
$graph->add_schema(namespace => 'QCGrind', label => 'Study',   unique => [qw(name)]);
my $john = $graph->add_node(namespace => 'QCGrind', label => 'Individual', properties => { name => 'John' });
my $jane = $graph->add_node(namespace => 'QCGrind', label => 'Individual', properties => { name => 'Jane' });
ok $graph->relate($jane, $sanger3, 'sampled'), 'relate() worked';
$graph->relate($john, $sanger2, 'sampled');
$graph->relate($john, $sanger1, 'sampled');
my $study = $graph->add_node(namespace => 'QCGrind', label => 'Study', properties => { name => 'Study of Disease_xyz' });
$graph->relate($study, $john, 'has_participant');
$graph->relate($study, $jane, 'has_participant');
my @samples = ($sanger1, $sanger2, $sanger3);

foreach my $i (1 .. 3) {
    my $library = $graph->add_node(namespace => 'QCGrind', label => 'Library', properties => { name => "Library$i" });
    $graph->relate($samples[$i - 1], $library, 'prepared');
    my $lane = $graph->add_node(namespace => 'QCGrind', label => 'Lane', properties => { name => "Lane$i", $i == 2 ? (foo => 'bar') : () });
    $graph->relate($library, $lane, 'sequenced');
}

@nodes = $graph->related_nodes($john, namespace => 'QCGrind', label => 'Sample');
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger1', 'sanger2'], 'related_nodes returned all Sample nodes related to John';
@nodes = $graph->related_nodes($john, namespace => 'QCGrind', label => 'Sample', properties => { sanger_id => 'sanger2' });
is_deeply [sort map { $_->get_property('sanger_id') } @nodes], ['sanger2'], 'related_nodes can be limited by params';
@nodes = $graph->related_nodes($john, namespace => 'QCGrind', label => 'Lane');
is_deeply [sort map { $_->get_property('name') } @nodes], [], 'related_nodes returned nothing looking for Lanes related to John because the default depth is only 1';
@nodes = $graph->related_nodes($john, namespace => 'QCGrind', label => 'Lane', max_depth => 3);
is_deeply [sort map { $_->get_property('name') } @nodes], [qw(Lane1 Lane2)], 'related_nodes returned the lanes related to John given a depth of 3';
@nodes = $graph->related_nodes($john, namespace => 'QCGrind', label => 'Lane', max_depth => 3, direction => '>');
is_deeply [sort map { $_->get_property('name') } @nodes], [qw(Lane1 Lane2)], 'related_nodes still works given an outgoing direction';
@nodes = $graph->related_nodes($john, namespace => 'QCGrind', label => 'Lane', max_depth => 3, direction => '<');
is_deeply [sort map { $_->get_property('name') } @nodes], [], 'related_nodes returns no lanes given an incoming direction';
@nodes = $graph->related_nodes($john, namespace => 'QCGrind', label => 'Study', direction => '<');
is_deeply [sort map { $_->get_property('name') } @nodes], ['Study of Disease_xyz'], 'related_nodes works with an incoming direction to find the study';
@nodes = $graph->related_nodes($john, namespace => 'QCGrind', label => 'Lane', max_depth => 3, direction => '>', properties => { foo => 'bar' });
is_deeply [sort map { $_->get_property('name') } @nodes], [qw(Lane2)], 'related_nodes still works when using max_depth, direction and properties all at the same time';

exit;
