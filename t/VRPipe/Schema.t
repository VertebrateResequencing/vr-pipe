#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 30;
    use VRPipeTest;
    use_ok('VRPipe::Schema');
}

ok my $schema = VRPipe::Schema->create('VRTrack'), 'able to create a schema instance';
ok my $graph = $schema->graph(), 'graph() method worked';
isa_ok($graph, 'VRPipe::Persistent::Graph');

ok my $sample = $schema->add('Sample', { name => 's1' }), 'add() method worked';
my $orig_sample_id = $sample->node_id();
$sample = $schema->add('Sample', { name => 's1', public_name => 'public_s1' });
is_deeply [$sample->{id}, $sample->{properties}], [$orig_sample_id, { name => 's1', public_name => 'public_s1' }], 'add() twice on the same unique property always does an update';
ok my @libs = $schema->add('Library', [{ id => 'l1' }, { id => 'l2' }], incoming => { type => 'prepared', node => $sample }), 'add() worked with incoming option, and for adding more than 1 at a time';

throws_ok { $schema->add('Foo', { foo => 'bar' }) } qr/'Foo' isn't a valid label for schema VRTrack/, 'add() throws when given an invalid label';
throws_ok { $schema->add('Sample', { id => '1' }) } qr/Parameter 'name' must be supplied/, 'add() throws when not given a required parameter';
throws_ok { $schema->add('Sample', { name => 's2', foo => 'bar' }) } qr/Property 'foo' supplied, but that isn't defined in the schema for VRTrack::Sample/, 'add() throws when given an invalid parameter';

ok my $lib1 = $schema->get('Library', { id => 'l1' }), 'get() method worked';
ok @libs = $schema->get('Library'), 'get() method worked with no properties arg';

$schema->add('Library', { id => 'l3' });
my $lib3 = $schema->get('Library', { id => 'l3' });
is $lib3->{properties}->{id}, 'l3', 'really added a library to the database';
$schema->delete($lib3);
$lib3 = $schema->get('Library', { id => 'l3' });
is $lib3, undef, 'delete() worked';

isa_ok($lib1, 'VRPipe::Schema::VRTrack::Library');
is $lib1->id, 'l1', 'auto-generated property method worked to get';
$lib1->name('libone');
is $lib1->name, 'libone', 'auto-generated property method worked to set';
$lib1 = $schema->get('Library', { id => 'l1' });
is $lib1->name, 'libone', 'the set was really in the database';

is_deeply $lib1->properties(flatten_parents => 1), { id => 'l1', name => 'libone', sample_name => 's1', sample_public_name => 'public_s1' }, 'properties() method worked with flatten_parents';
$lib1->add_properties({ name => 'lib1', tag => 'ATG' });
$lib1 = $schema->get('Library', { id => 'l1' });
is_deeply $lib1->properties(), { id => 'l1', name => 'lib1', tag => 'ATG' }, 'add_properties() method worked';
is $lib1->parent_property('sample_name'), 's1', 'parent_property() worked';

# test that we can get the latest data if another process updates a property
is $lib1->tag, 'ATG', 'tag starts out as ATG';
my $pid = fork;
if (defined $pid && $pid == 0) {
    $lib1->tag('ATGC');
    exit;
}
waitpid($pid, 0);
is $lib1->tag, 'ATG', 'tag still seems to be ATG after another process changed it';
$lib1->update_from_db;
is $lib1->tag, 'ATGC', 'tag is correct after using update_from_db()';

my @related = $lib1->related();
is_deeply [sort map { $_->node_id } @related], [$sample->node_id], 'related() worked';
$lib3 = $schema->add('Library', { id => 'l3' });
@related = $lib3->related();
is_deeply [sort map { $_->node_id } @related], [], 'related() returned nothing for a node with no relations';
$sample->relate_to($lib3, 'prepared');
@related = $lib3->related();
is_deeply [sort map { $_->node_id } @related], [$sample->node_id], 'relate_to() worked';
is $related[0]->name, 's1', 'related() returns working objects';
my $lane1 = $schema->add('Lane', { name => 'lane1' });
$lib3->relate_to($lane1, 'sequenced');
@related = $lane1->related(incoming => {});
is_deeply [sort map { $_->node_id } @related], [$lib3->node_id], 'related() worked with incoming specified';
$lib1->relate_to($lane1, 'sequenced', selfish => 1);
@related = $lane1->related(incoming => {});
is_deeply [sort map { $_->node_id } @related], [$lib1->node_id], 'relate_to(selfish => 1) worked';

exit;
