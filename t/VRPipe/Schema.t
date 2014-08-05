#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 52;
    use VRPipeTest;
    use_ok('VRPipe::Schema');
}

ok my $schema = VRPipe::Schema->create('VRTrack'), 'able to create a schema instance';
ok my $graph = $schema->graph(), 'graph() method worked';
isa_ok($graph, 'VRPipe::Persistent::Graph');

ok my $sample = $schema->add('Sample', { name => 's1' }), 'add() method worked';
my $orig_sample_id = $sample->node_id();
$sample = $schema->add('Sample', { name => 's1', public_name => 'public_sone' });
is_deeply [$sample->{id}, $sample->{properties}], [$orig_sample_id, { name => 's1', public_name => 'public_sone' }], 'add() twice on the same unique property does an update with labels that keep their history';
$sample = $schema->add('Sample', { name => 's1', public_name => 'public_s1' });
is_deeply [$sample->{id}, $sample->{properties}, $sample->changed()], [$orig_sample_id, { name => 's1', public_name => 'public_s1' }, { public_name => ['public_sone', 'public_s1'] }], 'add() again on the same unique property with a changed property does really update it, and we can find out what the previous value was';

# Sample is marked to store history but EBI_Submission is not, so repeat some of
# the above tests to cover both code paths
my $ebisub = $schema->add('EBI_Submission', { acc => 'ebi1', sub_date => 12345 });
my $orig_ebi_id = $ebisub->node_id();
$ebisub = $schema->add('EBI_Submission', { acc => 'ebi1', sub_date => 67890 });
is_deeply [$ebisub->{id}, $ebisub->{properties}, $ebisub->changed()], [$orig_ebi_id, { acc => 'ebi1', sub_date => 67890 }], 'add() twice on the same unique property with different other properties does an update on a history-less label';

ok my @libs = $schema->add('Library', [{ id => 'l1' }, { id => 'l2' }], incoming => { type => 'prepared', node => $sample }), 'add() worked with incoming option, and for adding more than 1 at a time';

throws_ok { $schema->add('Foo', { foo => 'bar' }) } qr/'Foo' isn't a valid label for schema VRTrack/, 'add() throws when given an invalid label';
throws_ok { $schema->add('Sample', { id => '1' }) } qr/Parameter 'name' must be supplied/, 'add() throws when not given a required parameter';
throws_ok { $schema->add('Sample', { name => 's2', foo => 'bar' }) } qr/Property 'foo' supplied, but that isn't defined in the schema for VRTrack::Sample/, 'add() throws when given an invalid parameter';
ok my $file = $schema->add('File', { path => '/path', foo => 'bar' }), 'arbitrary parameters can be supplied to a label defined with allow_anything => 1';
ok $file->add_properties({ cat => 'banana' }), 'add_properties() also worked with an arbitrary parameter';
is_deeply $file->{properties}, { path => '/path', foo => 'bar', cat => 'banana' }, 'We really do store whatever on a allow_anything label';

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

throws_ok { $sample->add_properties({ foo  => 'bar' }) } qr/Property 'foo' supplied, but that isn't defined in the schema for VRTrack::Sample/,           'add_properties() throws when given an invalid property';
throws_ok { $sample->add_properties({ name => 'bar' }) } qr/Property 'name' supplied, but that's unique for schema VRTrack::Sample and can't be changed/, 'add_properties() throws when given a unique property';
$sample->add_properties({ public_name => 'pn1', supplier_name => 'sn1', accession => 'acc1' });
is_deeply [$sample->properties(), $sample->changed()], [{ name => 's1', public_name => 'pn1', supplier_name => 'sn1', accession => 'acc1' }, { public_name => ['public_s1', 'pn1'], supplier_name => [undef, 'sn1'], accession => [undef, 'acc1'] }], 'add_properties() method worked and let us see what changed';
$sample->add_properties({ public_name => 'pn1', supplier_name => 'sup1' }, replace => 1);
is_deeply [$sample->properties(), $sample->changed()], [{ name => 's1', public_name => 'pn1', supplier_name => 'sup1' }, { supplier_name => ['sn1', 'sup1'], accession => ['acc1', undef] }], 'add_properties(replace => 1) method worked and let us see what changed';
$sample->add_properties({ public_name => 'pn1', supplier_name => 'sup1' });
is_deeply [$sample->properties(), $sample->changed()], [{ name => 's1', public_name => 'pn1', supplier_name => 'sup1' }], 'add_properties() with only old properties changes nothing';
throws_ok { $sample->name('bar') } qr/Property 'name' is unique for schema VRTrack::Sample and can't be changed/, 'property method name() throws when given a value to set';
is $sample->name(), 's1', 'but name() to get works fine';
$sample->supplier_name('supn1');
is_deeply [$sample->supplier_name(), $sample->properties(), $sample->changed()], ['supn1', { name => 's1', public_name => 'pn1', supplier_name => 'supn1' }, { supplier_name => ['sup1', 'supn1'] }], 'setting a property with the dynamically created method works and can tell us what we changed';

# test history
$sample = $schema->get('Sample', { name => 's1' });
my @history        = $sample->property_history();
my $timestamps_ok  = 0;
my $prev_timestamp = time();
my @history_properties;
foreach my $hist (@history) {
    my $timestamp = $hist->{timestamp};
    if ($timestamp <= $prev_timestamp) {
        $timestamps_ok++;
    }
    $prev_timestamp = $timestamp;
    
    push(@history_properties, $hist->{properties});
}
is_deeply [$timestamps_ok, @history_properties], [5, { public_name => 'pn1', supplier_name => 'supn1' }, { public_name => 'pn1', supplier_name => 'sup1' }, { public_name => 'pn1', supplier_name => 'sn1', accession => 'acc1' }, { public_name => 'public_s1' }, { public_name => 'public_sone' }], 'sample_property_history() gives us the complete history of properties over time';

# test its ok to set null values
$sample = $schema->add('Sample', { name => 's1', created_date => undef });
is_deeply [$sample->{id}, $sample->{properties}], [$orig_sample_id, { name => 's1', public_name => 'pn1', supplier_name => 'supn1' }], 'add() with null values for properties does not add those properties';
$sample = $schema->add('Sample', { name => 's1', created_date => 12 });
$sample = $schema->add('Sample', { name => 's1', created_date => undef });
is_deeply [$sample->{id}, $sample->{properties}], [$orig_sample_id, { name => 's1', public_name => 'pn1', supplier_name => 'supn1', created_date => 12 }], 'add() with null values for properties does not unset previously set properties';

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
my $lane1 = $schema->add('Lane', { unique => 'lane1', lane => 1 });
$lib3->relate_to($lane1, 'sequenced');
@related = $lane1->related(incoming => {});
is_deeply [sort map { $_->node_id } @related], [$lib3->node_id], 'related() worked with incoming specified';
$lib1->relate_to($lane1, 'sequenced', selfish => 1);
@related = $lane1->related(incoming => {});
is_deeply [sort map { $_->node_id } @related], [$lib1->node_id], 'relate_to(selfish => 1) worked';

# make sure we delete history nodes when we delete a schema node
$lib3->tag('A');
$lib3->tag('ATGC');
@related = $graph->related_nodes($lib3, outgoing => { max_depth => 500 }, return_history_nodes => 1);
is scalar(@related), 4, 'lib3 has 4 history nodes';
$schema->delete($lib3);
my $still_exist = 0;
foreach my $node (@related) {
    $still_exist++ if $graph->get_node_by_id($node->{id});
}
is $still_exist, 1, 'after deleting lib3, all the history nodes were also deleted, except for one used by another node';

# test Graph pass-through methods create_uuid() and date_to_epoch()
my $uuid = $schema->create_uuid();
like $uuid, qr/\w{8}-\w{4}-\w{4}-\w{4}-\w{12}/, 'create_uuid() worked';
is $schema->date_to_epoch('2013-05-10 06:45:32'), 1368168332, 'date_to_epoch() worked';

# unique uuid properties auto-fill if not supplied
ok my $bam_stats = $schema->add('Bam_Stats', { reads => 1, bases => 75 }), 'could add a new node without supplying its unique value when the unique is a uuid';
like $bam_stats->uuid, qr/\w{8}-\w{4}-\w{4}-\w{4}-\w{12}/, 'the resulting node has a uuid';

exit;
