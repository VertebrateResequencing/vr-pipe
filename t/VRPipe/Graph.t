#!/usr/bin/env perl
use strict;
use warnings;
use Parallel::ForkManager;
use Path::Class;

BEGIN {
    use Test::Most tests => 42;
    use VRPipeTest;
    use_ok('VRPipe::Persistent::Graph');
}

ok my $graph = VRPipe::Persistent::Graph->new, 'able to create a Graph object';
isa_ok($graph, 'VRPipe::Persistent::Graph');

ok $graph->drop_database, 'we were able to drop the test database';

ok $graph->add_schema(namespace => 'VRTrack', label => 'Sample', unique => [qw(sanger_id)], indexed => [qw(uuid public_name)]), 'add_schema worked';
is $graph->add_schema(namespace => 'VRTrack', label => 'Sample', unique => [qw(sanger_id)], indexed => [qw(uuid public_name)]), 0, 'add_schema a second time with same args does nothing';
throws_ok { $graph->add_schema(namespace => 'OtherNS', label => 'Sam|ple', unique => [qw(name)]) } qr/neither namespace or label may contain the \| character/, 'add_schema throws if you try a label with the pipe character in it';
throws_ok { $graph->add_schema(namespace => 'OtherNS', label => 'Sample', unique => ['nam|e']) } qr/FATAL ERROR.+parameter may not contain the \| character/sm, 'add_schema throws if you try a parameter with the pipe character in it';
$graph->throw_with_no_stacktrace(1);
throws_ok { $graph->add_schema(namespace => 'OtherNS', label => 'Sample', unique => ['nam|e']) } qr/^parameter may not contain the \| character/, 'throw has no stacktrace if throw_with_no_stacktrace is turned on';
ok $graph->add_schema(namespace => 'OtherNS', label => 'Sample', unique => [qw(name)]), 'add_schema worked for the same label in a different namespace';
my ($unique, $indexed) = $graph->get_schema(namespace => 'VRTrack', label => 'Sample');
is_deeply [$unique, $indexed], [[qw(sanger_id)], [qw(uuid public_name)]], 'get_schema returned the unique and indexed fields from its internal cache';
$VRPipe::Persistent::Graph::schemas = {};
($unique, $indexed) = $graph->get_schema(namespace => 'VRTrack', label => 'Sample');
is_deeply [$unique, $indexed], [[qw(sanger_id)], [qw(uuid public_name)]], 'get_schema returned the unique and indexed fields from the database';

ok my $node = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' }), 'was able to add a node';
my $orig_node_id = $graph->node_id($node);
my $sanger1 = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' });
is $graph->node_id($sanger1), $orig_node_id, 'calling add_node twice with the same args returns the same node';
throws_ok { $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', public_name => 'public2' }) } qr/already exists with label vdt.+\|VRTrack\|Sample/, 'add_node throws when used twice with the same unique arg but different other arg';
ok $graph->delete_node($sanger1), 'delete_node() worked';
undef $sanger1;
($sanger1) = $graph->get_nodes(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1' });
ok !$sanger1, 'the node really is gone';
$sanger1 = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' });
isnt $graph->node_id($sanger1), $orig_node_id, 'calling add_node again with the same args after deleting returns a new node';
throws_ok { $graph->add_node(namespace => 'VRTrack', label => 'Sampl', properties => { sanger_id => 'sanger1', uuid => 'uuuuu', public_name => 'public1' }) } qr/You must first create a schema for namespace/, 'add_node throws when supplied with a label for which there is no schema';

# get nodes by their properties (first add some more for testing purposes)
my $sanger2 = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger2', public_name => 'public1' });
my $sanger3 = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger3', public_name => 'public2' });
my @nodes = $graph->get_nodes(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1' });
is_deeply [sort map { $graph->node_property($_, 'sanger_id') } @nodes], ['sanger1'], 'get_nodes returned expected node when searching on a unique value';
@nodes = $graph->get_nodes(namespace => 'VRTrack', label => 'Sample', properties => { public_name => 'public1' });
is_deeply [sort map { $graph->node_property($_, 'sanger_id') } @nodes], ['sanger1', 'sanger2'], 'get_nodes returned expected nodes when searching on a non-unique value';

# relationship tests (first add even more nodes to test with)
ok $graph->add_schema(namespace => 'VRTrack', label => 'Individual', unique => [qw(name)]), 'add_schema worked with different args and no indexed';
$graph->add_schema(namespace => 'VRTrack', label => 'Library', unique => [qw(name)]);
$graph->add_schema(namespace => 'VRTrack', label => 'Lane',    unique => [qw(name)]);
$graph->add_schema(namespace => 'VRTrack', label => 'Study',   unique => [qw(name)]);
my $john = $graph->add_node(namespace => 'VRTrack', label => 'Individual', properties => { name => 'John' });
my $jane = $graph->add_node(namespace => 'VRTrack', label => 'Individual', properties => { name => 'Jane' });
ok $graph->relate($jane, $sanger3, type => 'sampled'), 'relate() worked';
$graph->relate($john, $sanger2, type => 'sampled');
$graph->relate($john, $sanger1, type => 'sampled');
my $study = $graph->add_node(namespace => 'VRTrack', label => 'Study', properties => { name => 'Study of Disease_xyz' });
my @root_ids = ($graph->node_id($study));
$graph->relate($study, $john, type => 'has_participant');
$graph->relate($study, $jane, type => 'has_participant');
my @samples = ($sanger1, $sanger2, $sanger3);

foreach my $i (1 .. 3) {
    my $library = $graph->add_node(namespace => 'VRTrack', label => 'Library', properties => { name => "Library$i" });
    $graph->relate($samples[$i - 1], $library, type => 'prepared');
    my $lane = $graph->add_node(namespace => 'VRTrack', label => 'Lane', properties => { name => "Lane$i", $i == 2 ? (foo => 'bar') : () });
    $graph->relate($library, $lane, type => 'sequenced');
}

@nodes = $graph->related_nodes($john, undirected => { namespace => 'VRTrack', label => 'Sample' });
is_deeply [sort map { $graph->node_property($_, 'sanger_id') } @nodes], ['sanger1', 'sanger2'], 'related_nodes returned all Sample nodes related to John';
@nodes = $graph->related_nodes($john, undirected => { namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger2' } });
is_deeply [sort map { $graph->node_property($_, 'sanger_id') } @nodes], ['sanger2'], 'related_nodes can be limited by params';
@nodes = $graph->related_nodes($john, undirected => { namespace => 'VRTrack', label => 'Lane' });
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], [], 'related_nodes returned nothing looking for Lanes related to John because the default depth is only 1';
@nodes = $graph->related_nodes($john, undirected => { namespace => 'VRTrack', label => 'Lane', max_depth => 3 });
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], [qw(Lane1 Lane2)], 'related_nodes returned the lanes related to John given a depth of 3';
@nodes = $graph->related_nodes($john, outgoing => { namespace => 'VRTrack', label => 'Lane', max_depth => 3 });
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], [qw(Lane1 Lane2)], 'related_nodes still works given an outgoing direction';
@nodes = $graph->related_nodes($john, incoming => { namespace => 'VRTrack', label => 'Lane', max_depth => 3 });
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], [], 'related_nodes returns no lanes given an incoming direction';
@nodes = $graph->related_nodes($john, incoming => { namespace => 'VRTrack', label => 'Study' });
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], ['Study of Disease_xyz'], 'related_nodes works with an incoming direction to find the study';
@nodes = $graph->related_nodes($john, outgoing => { namespace => 'VRTrack', label => 'Lane', max_depth => 3, properties => { foo => 'bar' } });
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], [qw(Lane2)], 'related_nodes still works when using max_depth, direction and properties all at the same time';
@nodes = $graph->related_nodes($john, incoming => {}, outgoing => { min_depth => 3 });
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], ['Lane1', 'Lane2', 'Study of Disease_xyz'], 'related_nodes works with both directions and min_depth';
@nodes = $graph->related_nodes($sanger2);
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], [qw(John Library2)], 'related_nodes defaults to giving all nodes 1 step away';

# add some more nodes to test the visualisation and add_nodes() bulk creation
# with relationships at the same time
$study = $graph->add_node(namespace => 'VRTrack', label => 'Study', properties => { name => 'Study of Disease_abc' });
push(@root_ids, $graph->node_id($study));
$graph->relate($study, $jane, type => 'has_participant');
$graph->add_schema(namespace => 'OtherNS', label => 'Workplace', unique => [qw(name)], indexed => []);
my $workplace = $graph->add_node(namespace => 'OtherNS', label => 'Workplace', properties => { name => 'Sanger' });
$graph->add_schema(namespace => 'OtherNS', label => 'Individual', unique => [qw(name)], indexed => []);
my @props;

for (1 .. 1000) {
    push(@props, { name => 'Person' . $_ });
}
@nodes = $graph->add_nodes(namespace => 'OtherNS', label => 'Individual', properties => \@props, incoming => { node => $workplace, type => 'has_employee' }); # ~9 seconds; is that good or bad?
is scalar(@nodes), 1000, 'add_nodes() worked and returned 1000 nodes';
@props = ();
for (951 .. 1050) {
    push(@props, { name => 'Person' . $_ });
}
@nodes = $graph->add_nodes(namespace => 'OtherNS', label => 'Individual', properties => \@props, incoming => { node => $workplace, type => 'has_employee' });
is scalar(@nodes), 100, 'add_nodes() worked and returned 100 nodes when we requested to create some already existing nodes';
my $employee = $nodes[0];
@nodes = $graph->related_nodes($workplace);
is scalar(@nodes), 1050, 'add_nodes() incoming option worked and related all 1050 nodes to another node';
my $campus = $graph->add_node(namespace => 'OtherNS', label => 'Workplace', properties => { name => 'Genome Campus' }, outgoing => { node => $workplace, type => 'contains' });
push(@root_ids, $graph->node_id($campus));
($node) = $graph->related_nodes($workplace, incoming => {});
is $graph->node_property($node, 'name'), 'Genome Campus', 'add_node() outgoing option worked and created the relationship';
my $ebi = $graph->add_node(namespace => 'OtherNS', label => 'Workplace', properties => { name => 'EBI' }, incoming => { node => $campus, type => 'contains' }, outgoing => { node => $employee, type => 'has_employee' });
@nodes = $graph->related_nodes($ebi);
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], ['Genome Campus', $graph->node_property($employee, 'name')], 'add_node() worked with both incoming and outgoing specified at once';

# add nodes to test display of images, and test the uuid creation helper method
# and required schema option
my ($lane1) = $graph->get_nodes(namespace => 'VRTrack', label => 'Lane', properties => { name => 'Lane1' });
$graph->add_schema(namespace => 'VRPipe', label => 'StepResult', unique => [qw(uuid)]);
$graph->add_schema(namespace => 'VRPipe', label => 'Image', unique => [qw(path)], required => ['type']);
my $step_result = $graph->add_node(namespace => 'VRPipe', label => 'StepResult', properties => { uuid => $graph->create_uuid });
push(@root_ids, $graph->node_id($step_result));
like $graph->node_property($step_result, 'uuid'), qr/\w{8}-\w{4}-\w{4}-\w{4}-\w{12}/, 'create_uuid() method worked to populate a property';
throws_ok { $graph->add_node(namespace => 'VRPipe', label => 'StepResult', properties => { foo => 'bar' }) } qr/Parameter 'uuid' must be supplied/, 'add_node throws when a unique parameter is not supplied';
my $image = $graph->add_node(namespace => 'VRPipe', label => 'Image', properties => { path => file(qw(t data qcgraph.png))->absolute->stringify, type => 'png' }, incoming => { node => $step_result, type => 'has_file' });
throws_ok { $graph->add_node(namespace => 'VRPipe', label => 'Image', properties => { path => '/foo' }) } qr/Parameter 'type' must be supplied/, 'add_node throws when a specified required parameter is not supplied';

# test root_nodes method
@nodes = $graph->root_nodes();
is_deeply {
    map { $graph->node_id($_) => 1 } @nodes;
}, { map { $_ => 1 } @root_ids }, 'root_nodes() worked as execpted';

exit;
