#!/usr/bin/env perl
use strict;
use warnings;
use Parallel::ForkManager;
use Path::Class;

BEGIN {
    use Test::Most tests => 84;
    use VRPipeTest;
    use_ok('VRPipe::Persistent::Graph');
}

ok my $graph = VRPipe::Persistent::Graph->new, 'able to create a Graph object';
isa_ok($graph, 'VRPipe::Persistent::Graph');

ok $graph->drop_database, 'we were able to drop the test database';

ok $graph->add_schema(namespace => 'VRTrack', label => 'Sample', unique => [qw(sanger_id)], indexed => [qw(uuid public_name)]),               'add_schema worked';
is $graph->add_schema(namespace => 'VRTrack', label => 'Sample', unique => [qw(sanger_id)], indexed => [qw(uuid public_name)]),               0, 'add_schema a second time with same args does nothing';
ok $graph->add_schema(namespace => 'VRTrack', label => 'Sample', unique => [qw(sanger_id)], indexed => [qw(uuid public_name supplier_name)]), 'add_schema updated an existing schema';
throws_ok { $graph->add_schema(namespace => 'OtherNS', label => 'Sam|ple', unique => [qw(name)]) } qr/neither namespace or label may contain the \| character/, 'add_schema throws if you try a label with the pipe character in it';
throws_ok { $graph->add_schema(namespace => 'OtherNS', label => 'Sample', unique => ['nam|e']) } qr/FATAL ERROR.+parameter may not contain the \| character/sm, 'add_schema throws if you try a parameter with the pipe character in it';
$graph->throw_with_no_stacktrace(1);
throws_ok { $graph->add_schema(namespace => 'OtherNS', label => 'Sample', unique => ['nam|e']) } qr/^parameter may not contain the \| character/, 'throw has no stacktrace if throw_with_no_stacktrace is turned on';
ok $graph->add_schema(namespace => 'OtherNS', label => 'Sample', unique => [qw(name)]), 'add_schema worked for the same label in a different namespace';
my ($unique, $indexed) = $graph->get_schema(namespace => 'VRTrack', label => 'Sample');
is_deeply [$unique, $indexed], [[qw(sanger_id)], [qw(uuid public_name supplier_name)]], 'get_schema returned the unique and indexed fields from its internal cache';
$VRPipe::Persistent::Graph::schemas = {};
($unique, $indexed) = $graph->get_schema(namespace => 'VRTrack', label => 'Sample');
is_deeply [$unique, $indexed], [[qw(sanger_id)], [qw(uuid public_name supplier_name)]], 'get_schema returned the unique and indexed fields from the database';

ok my $node = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', public_name => 'publicOne' }), 'was able to add a node';
my $orig_node_id = $graph->node_id($node);
my $sanger1 = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', public_name => 'publicOne' });
is_deeply [$sanger1->{id}, $sanger1->{properties}], [$orig_node_id, { sanger_id => 'sanger1', public_name => 'publicOne' }], 'calling add_node twice with the same args returns the same node, which has expected properties';
throws_ok { $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', public_name => 'public2' }) } qr/already exists with label vdt.+\|VRTrack\|Sample/, 'add_node throws when used twice with the same unique arg but different other arg';
ok $sanger1 = $graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'sanger1', public_name => 'public1', uuid => 'uuuuu' }, update => 1), 'add_node(update => 1) works when used twice with the same unique arg and different other args';
is_deeply [$sanger1->{id}, $sanger1->{properties}], [$orig_node_id, { sanger_id => 'sanger1', public_name => 'public1', uuid => 'uuuuu' }], 'calling add_node in update mode twice with the different non-unique args returns the same node, which has updated properties';
ok $graph->delete_node($sanger1), 'delete_node() worked';
$VRPipe::Persistent::Graph::schemas = {};
($unique, $indexed) = $graph->get_schema(namespace => 'VRTrack', label => 'Sample');
ok $unique, 'delete_node() did not delete all nodes indiscriminately';
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
@nodes = $graph->related_nodes($study, incoming => { min_depth => 0 }, outgoing => { min_depth => 0 });
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], [qw(Jane John)], 'related_nodes works correctly with a min_depth of 0';

# test deleting relationships
my $samson  = $graph->add_node(namespace => 'VRTrack', label => 'Individual', properties => { name => 'Samson' });
my $delilah = $graph->add_node(namespace => 'VRTrack', label => 'Individual', properties => { name => 'Delilah' });
$graph->relate($samson,  $delilah, type => 'loves');
$graph->relate($delilah, $samson,  type => 'betrays');
@nodes = $graph->related_nodes($samson, outgoing => {});
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], [qw(Delilah)], 'Samson is related to Delilah';
@nodes = $graph->related_nodes($delilah, outgoing => {});
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], [qw(Samson)], 'Delilah is related to Samson';
$graph->divorce($samson, $delilah, type => 'loves');
@nodes = $graph->related_nodes($samson, outgoing => {});
is scalar(@nodes), 0, 'after divorcing on loves, Samson no longer related to Delilah in that direction';
@nodes = $graph->related_nodes($delilah, outgoing => {});
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], [qw(Samson)], 'Delilah remains related to Samson';
$graph->divorce($samson, $delilah);
@nodes = $graph->related_nodes($samson, outgoing => {});
is scalar(@nodes), 0, 'after divorcing without specifying type, Samson remains unrelated to Delilah in that direction';
@nodes = $graph->related_nodes($delilah, outgoing => {});
is scalar(@nodes), 0, 'and now Delilah is not related to Samson';
$graph->delete_node($samson);
$graph->delete_node($delilah);

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
my ($tom, $jim) = $graph->add_nodes(namespace => 'OtherNS', label => 'Individual', properties => [{ name => 'Tom' }, { name => 'Jim' }], incoming => [{ type => 'has_employee', node_spec => { namespace => 'OtherNS', label => 'Workplace', properties => { name => 'Genome Campus' } } }, { type => 'has_employee', node_spec => { namespace => 'OtherNS', label => 'Workplace', properties => { name => 'EBI' } } }]);
@nodes = $graph->related_nodes($jim);
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], ['EBI', 'Genome Campus'], 'add_nodes() worked with 2 incoming nodes defined with specs';
@nodes = $graph->related_nodes($ebi, outgoing => {});
is scalar(@nodes), 3, 'before creating mass relationships, had just 3';
my @rel_details;

for (1 .. 1000) {
    push(@rel_details, { from => $ebi, to => { namespace => 'OtherNS', label => 'Individual', properties => { name => 'Person' . $_ } }, type => 'has_employee' });
}
$graph->create_mass_relationships(\@rel_details);
@nodes = $graph->related_nodes($ebi, outgoing => {});
is scalar(@nodes), 1002, 'after creating mass relationships, had 1002';

# add nodes to test display of images, and test the uuid creation helper method
# and required schema option
my ($lane1) = $graph->get_nodes(namespace => 'VRTrack', label => 'Lane', properties => { name => 'Lane1' });
$graph->add_schema(namespace => 'VRPipe', label => 'StepResult', unique => [qw(uuid)]);
$graph->add_schema(namespace => 'VRPipe', label => 'Image', unique => [qw(path)], required => ['type']);
my $uuid = $graph->create_uuid;
my $step_result = $graph->add_node(namespace => 'VRPipe', label => 'StepResult', properties => { uuid => $uuid }, incoming => { node => $lane1, type => 'has_result' });
like $graph->node_property($step_result, 'uuid'), qr/\w{8}-\w{4}-\w{4}-\w{4}-\w{12}/, 'create_uuid() method worked to populate a property';
throws_ok { $graph->add_node(namespace => 'VRPipe', label => 'StepResult', properties => { foo => 'bar' }) } qr/Parameter 'uuid' must be supplied/, 'add_node throws when a unique parameter is not supplied';
my $image = $graph->add_node(namespace => 'VRPipe', label => 'Image', properties => { path => file(qw(t data qcgraph.png))->absolute->stringify, type => 'png' }, incoming => { node => $step_result, type => 'has_file' });
throws_ok { $graph->add_node(namespace => 'VRPipe', label => 'Image', properties => { path => '/foo' }) } qr/Parameter 'type' must be supplied/, 'add_node throws when a specified required parameter is not supplied';

# test root_nodes method
@nodes = $graph->root_nodes();
my %node_ids = map { $graph->node_id($_) => 1 } @nodes;
is_deeply \%node_ids, { map { $_ => 1 } @root_ids }, 'root_nodes() worked as execpted';

# make sure there are no problems getting a ton of nodes back from a cypher
# query
my $data = $graph->_run_cypher([['match n return n']]);
is scalar(@{ $data->{nodes} }), 1070, 'no problems returning lots of nodes';
my $encoded = $graph->json_encode($data);
my $decoded = $graph->json_decode($encoded);
is scalar(@{ $decoded->{nodes} }), 1070, 'json_encode/decode successfully roundtrips on lots of nodes';

# test the check_parents option to node_property
is $graph->node_property($image, 'vrtrack_lane_name'), undef, 'by default, parent properties are not available';
is_deeply $graph->node_properties($image), { path => file(qw(t data qcgraph.png))->absolute->stringify, type => 'png' }, 'node_properties just gives image details';
is $graph->node_property($image, 'vrtrack_lane_name', check_parents => 1), 'Lane1', 'in check_parents mode we can access a lane detail';
is_deeply $graph->node_properties($image), { path => file(qw(t data qcgraph.png))->absolute->stringify, type => 'png' }, 'node_properties still just gives image details';
is_deeply $graph->node_properties($image, flatten_parents => 1), { path => file(qw(t data qcgraph.png))->absolute->stringify, type => 'png', stepresult_uuid => $uuid, vrtrack_lane_name => 'Lane1', vrtrack_library_name => 'Library1', vrtrack_sample_sanger_id => 'sanger1', vrtrack_sample_public_name => 'public1', vrtrack_sample_uuid => 'uuuuu', vrtrack_individual_name => 'John', vrtrack_study_name => 'Study of Disease_xyz' }, 'in flatten_parents mode we see all parental properties';

# we can change existing properties and add new ones and remove them
$graph->node_add_properties($step_result, { foo => 'bar', cat => 'dog' });
is_deeply $step_result->{properties}, { uuid => $uuid, foo => 'bar', cat => 'dog' }, 'node_add_properties() adds properties correctly';
$graph->node_add_properties($step_result, { foo => 'baz', lemur => 'llama' });
my ($fresh_step_result) = $graph->get_nodes(namespace => 'VRPipe', label => 'StepResult', properties => { uuid => $uuid });
is_deeply $fresh_step_result->{properties}, { uuid => $uuid, foo => 'baz', cat => 'dog', lemur => 'llama' }, 'node_add_properties() adds and changes properties correctly, and the results are really in the database';
$graph->node_set_properties($fresh_step_result, { uuid => $uuid, foo => 'baz', lemur => 'lamella' });
is_deeply $fresh_step_result->{properties}, { uuid => $uuid, foo => 'baz', lemur => 'lamella' }, 'node_set_properties() updates properties correctly and removes unspecified ones';
($fresh_step_result) = $graph->get_nodes(namespace => 'VRPipe', label => 'StepResult', properties => { uuid => $uuid });
is_deeply $fresh_step_result->{properties}, { uuid => $uuid, foo => 'baz', lemur => 'lamella' }, 'node_set_properties() had its effect in the database';
$graph->node_add_properties($step_result, { bad_property => 'foo' });
($fresh_step_result) = $graph->get_nodes(namespace => 'VRPipe', label => 'StepResult', properties => { uuid => $uuid });
is_deeply $fresh_step_result->{properties}, { uuid => $uuid, foo => 'baz', lemur => 'lamella', bad_property => 'foo' }, 'added a bad property...';
$graph->node_remove_property($step_result, 'bad_property');
is_deeply $step_result->{properties}, { uuid => $uuid, foo => 'baz', lemur => 'lamella' }, 'node_remove_property() worked to remove that property';
($fresh_step_result) = $graph->get_nodes(namespace => 'VRPipe', label => 'StepResult', properties => { uuid => $uuid });
is_deeply $fresh_step_result->{properties}, { uuid => $uuid, foo => 'baz', lemur => 'lamella' }, 'and it is really removed from the database';

# we can get a node by its database id
$node = $graph->get_node_by_id($graph->node_id($step_result));
is_deeply $node->{properties}, { uuid => $uuid, foo => 'baz', lemur => 'lamella' }, 'get_node_by_id() worked';

# usually you can have node related to an unlimited number of other nodes,
# but sometimes you want it to only connect to a single node of a certain type,
# and an update should remove any existing links before applying the new one
my ($library1) = $graph->get_nodes(namespace => 'VRTrack', label => 'Library', properties => { name => "Library1" });
my ($library2) = $graph->get_nodes(namespace => 'VRTrack', label => 'Library', properties => { name => "Library2" });
$graph->relate($library2, $lane1, type => 'sequenced');
@nodes = $graph->related_nodes($lane1, incoming => { namespace => 'VRTrack', label => 'Library' });
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], ['Library1', 'Library2'], 'using standard relate(), a lane can belong to more than 1 incoming node';
$graph->relate($library1, $lane1, type => 'sequenced', selfish => 1);
@nodes = $graph->related_nodes($lane1, incoming => { namespace => 'VRTrack', label => 'Library' });
is_deeply [sort map { $graph->node_property($_, 'name') } @nodes], ['Library1'], 'using relate(selfish => 1), a lane only belongs to 1 incoming node';
my $image2 = $graph->add_node(namespace => 'VRPipe', label => 'Image', properties => { path => 'img2', type => 'png' }, incoming => { node => $image, type => 'sub_image' });
my $image3 = $graph->add_node(namespace => 'VRPipe', label => 'Image', properties => { path => 'img3', type => 'png' });
$graph->relate($image, $image3, type => 'sub_image');
@nodes = $graph->related_nodes($image, outgoing => {});
is_deeply [sort map { $graph->node_property($_, 'path') } @nodes], ['img2', 'img3'], 'using standard relate() we can have a node connected to 2 outgoing nodes';
$graph->relate($image, $image2, type => 'sub_image', replace => 1);
@nodes = $graph->related_nodes($image, outgoing => {});
is_deeply [sort map { $graph->node_property($_, 'path') } @nodes], ['img2'], 'using relate(replace => 1) we end up with only 1 outgoing node';
$graph->delete_node($image2);
$graph->delete_node($image3);

# make sure we can add_node(update) when we have required props
$graph->add_schema(namespace => 'Foo', label => 'Bar', unique => [qw(id)], required => [qw(name)]);
ok $sanger1 = $graph->add_node(namespace => 'Foo', label => 'Bar', properties => { id => 'f1', name => 'name1' }, update => 1), 'add_node(update => 1) works when used on a schema with required properties';
$graph->drop_schema(namespace => 'Foo', label => 'Bar');

# test the date_to_epoch utility method
is $graph->date_to_epoch('2013-05-10 06:45:32'), 1368168332, 'date_to_epoch() worked';

# test adding nodes in enqueue mode
ok !$graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'enqueue1' }, enqueue => 1), 'add_node(enqueue => 1) returns nothing';
ok !$graph->get_nodes(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'enqueue1' }), 'it did not add a node to the database';
ok !$graph->add_node(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'enqueue2' }, enqueue => 1), 'add_node(enqueue => 1) worked again';
ok my @queued = $graph->dispatch_queue(), 'could call dispatch_queue()';
is_deeply [sort map { $graph->node_property($_, 'sanger_id') } @queued], ['enqueue1', 'enqueue2'], 'dispatch_queue() returned the enqueued nodes';
is_deeply [sort map { $graph->node_property($_, 'sanger_id') } ($graph->get_nodes(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'enqueue1' }), $graph->get_nodes(namespace => 'VRTrack', label => 'Sample', properties => { sanger_id => 'enqueue2' }))], ['enqueue1', 'enqueue2'], 'after dispatch_queue() the enqueued nodes are in the database';

# we can set properties on a relationship (this is so far a rare-case usage, so
# there isn't an easy way to set this at node/relationship creation time, nor
# easy way to get relationships)
$image2 = $graph->add_node(namespace => 'VRPipe', label => 'Image', properties => { path => 'img2', type => 'png' }, incoming => { node => $image, type => 'sub_image' });
my ($rel) = @{ $graph->related($image2, undef, {})->{relationships} };
$graph->relationship_set_properties($rel, { reason => 'likes subs' });
($rel) = @{ $graph->related($image2, undef, {})->{relationships} };
is_deeply $rel->{properties}, { reason => 'likes subs' }, 'relationship_set_properties() worked';

exit;
