#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 169;
    use VRPipeTest;
    use_ok('VRPipe::Schema');
    use_ok('VRPipe::File');
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
is $sample->unique_property, 's1', 'unique_property() works';

# Sample is marked to store history but EBI_Submission is not, so repeat some of
# the above tests to cover both code paths
my $ebisub = $schema->add('EBI_Submission', { acc => 'ebi1', sub_date => 12345 });
my $orig_ebi_id = $ebisub->node_id();
$ebisub = $schema->add('EBI_Submission', { acc => 'ebi1', sub_date => 67890 });
is_deeply [$ebisub->{id}, $ebisub->{properties}, $ebisub->changed()], [$orig_ebi_id, { acc => 'ebi1', sub_date => 67890 }], 'add() twice on the same unique property with different other properties does an update on a history-less label';

ok my @libs = $schema->add('Library', [{ id => 'l1', center_name => 'sanger1' }, { id => 'l2' }], incoming => { type => 'prepared', node => $sample }), 'add() worked with incoming option, and for adding more than 1 at a time';

throws_ok { $schema->add('Foo', { foo => 'bar' }) } qr/'Foo' isn't a valid label for schema VRTrack/, 'add() throws when given an invalid label';
throws_ok { $schema->add('Sample', { id => '1' }) } qr/Parameter 'name' must be supplied/, 'add() throws when not given a required parameter';
throws_ok { $schema->add('Sample', { name => 's2', foo => 'bar' }) } qr/Property 'foo' supplied, but that isn't defined in the schema for VRTrack::Sample/, 'add() throws when given an invalid parameter';
my $bs_date = time();
ok my $bs = $schema->add('Bam_Stats', { uuid => 'uuid', mode => 'mode', options => 'opts', 'raw total sequences' => 100, foo => 'bar', date => $bs_date }), 'arbitrary parameters can be supplied to a label defined with allow_anything => 1';
ok $bs->add_properties({ cat => 'banana' }), 'add_properties() also worked with an arbitrary parameter';
is_deeply $bs->{properties}, { uuid => 'uuid', mode => 'mode', options => 'opts', 'raw total sequences' => 100, foo => 'bar', cat => 'banana', date => $bs_date }, 'We really do store whatever on a allow_anything label';

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

is_deeply $lib1->properties(flatten_parents => 1), { id => 'l1', name => 'libone', sample_name => 's1', sample_public_name => 'public_s1', center_name => 'sanger1' }, 'properties() method worked with flatten_parents';
$lib1->add_properties({ name => 'lib1', tag => 'ATG' });
$lib1 = $schema->get('Library', { id => 'l1' });
is_deeply $lib1->properties(), { id => 'l1', name => 'lib1', center_name => 'sanger1', tag => 'ATG', }, 'add_properties() method worked';
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
# also check it works on allow_anything labels
my $fooless_bs = $schema->add('Bam_Stats', { uuid => 'uuid2', mode => 'mode', options => 'opts', 'raw total sequences' => 100, foo => undef, date => $bs_date });
is_deeply $fooless_bs->{properties}, { uuid => 'uuid2', mode => 'mode', options => 'opts', 'raw total sequences' => 100, date => $bs_date }, 'add() with null values for properties does not add those properties for an allow_anything label';

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
my $lib4 = $schema->add('Library', { id => 'l4', center_name => 'sanger' });
$sample->relate_to($lib4, 'prepared');
my $lane1 = $schema->add('Lane', { unique => 'lane1', lane => 1 });
$lib3->relate_to($lane1, 'sequenced');
@related = $lane1->related(incoming => {});
is_deeply [sort map { $_->node_id } @related], [$lib3->node_id], 'related() worked with incoming specified';
$lib1->relate_to($lane1, 'sequenced', selfish => 1, properties => { machine => 'big_one' });
@related = $lane1->related(incoming => {});
is_deeply [sort map { $_->node_id } @related], [$lib1->node_id], 'relate_to(selfish => 1) worked';
my ($r) = @{ $graph->related($lib1, undef, undef, { type => 'sequenced' })->{relationships} };
is $r->{properties}->{machine}, 'big_one', 'relate_to(properties => {}) worked';

my $closest = $sample->closest('VRTrack', 'Library', direction => 'outgoing');
like $closest->id, qr/l[12]/, 'closest(outgoing) worked';
my @closest = $sample->closest('VRTrack', 'Library', direction => 'outgoing', all => 1);
is_deeply [sort map { $_->id } @closest], ['l1', 'l2', 'l3', 'l4'], 'closest(all) worked';
@closest = $sample->closest('VRTrack', 'Library', direction => 'outgoing', all => 1, properties => [['center_name', 'sanger']]);
is_deeply [sort map { $_->id } @closest], ['l4'], 'closest(all) with a specified property limit worked';
@closest = $sample->closest('VRTrack', 'Library', direction => 'outgoing', all => 1, properties => [['center_name', 'san\wer.*', 1]]);
is_deeply [sort map { $_->id } @closest], ['l1', 'l4'], 'closest(all) with a regex property limit worked';
@closest = $sample->closest('VRTrack', 'Library', direction => 'outgoing', all => 1, properties => [['center_name', 'san\wer.*', 1], ['id', 'l[12]', 1]]);
is_deeply [sort map { $_->id } @closest], ['l1'], 'closest(all) with 2 regex property limits worked';
$closest = $lib3->closest('VRTrack', 'Sample', direction => 'outgoing');
is $closest, undef, 'closest(outgoing) returns nothing when searching for non-existing';
$closest = $lib3->closest('VRTrack', 'Sample', direction => 'incoming');
is $closest->name, 's1', 'closest(incoming) worked';
$closest = $lib3->closest('VRTrack', 'Sample');
is $closest->name, 's1', 'closest() worked';

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
ok my $bam_stats = $schema->add('Bam_Stats', { mode => 'normal', options => '-foo', 'raw total sequences' => 100, date => $bs_date }), 'could add a new node without supplying its unique value when the unique is a uuid';
like $bam_stats->uuid, qr/\w{8}-\w{4}-\w{4}-\w{4}-\w{12}/, 'the resulting node has a uuid';
is $bam_stats->property('raw total sequences'), 100, 'we can get a property that has spaces in the name';

# test the VRTrack-specific ensure_sequencing_hierarchy method
ok my $hierarchy = $schema->ensure_sequencing_hierarchy(lane => 'esh_lane1', library => 'esh_library1', sample => 'esh_sample1', study => 'esh_study1', group => 'esh_group1', taxon => 'esh_taxon1'), 'ensure_sequencing_hierarchy() worked';
ok my $eshlane = $schema->get('Lane', { unique => 'esh_lane1' }), 'ensure_sequencing_hierarchy() created a Lane';
my $eshlib = $schema->get('Library', { id   => 'esh_library1' });
my $eshsam = $schema->get('Sample',  { name => 'esh_sample1' });
my $eshstu = $schema->get('Study',   { id   => 'esh_study1' });
my $eshgro = $schema->get('Group',   { name => 'esh_group1' });
my $eshtax = $schema->get('Taxon',   { id   => 'esh_taxon1' });
my %hierarchy_props = map { $_ => $hierarchy->{$_}->properties } keys %{$hierarchy};
is_deeply [\%hierarchy_props, [sort { $a <=> $b } map { $_->node_id } $eshlane->related(incoming => { max_depth => 10 })]], [{ lane => { unique => 'esh_lane1', lane => 'esh_lane1' }, library => { id => 'esh_library1' }, sample => { name => 'esh_sample1' }, study => { id => 'esh_study1' }, group => { name => 'esh_group1' }, taxon => { id => 'esh_taxon1' } }, [$eshlib->node_id, $eshsam->node_id, $eshstu->node_id, $eshgro->node_id, $eshtax->node_id]], 'ensure_sequencing_hierarchy() created the whole hierarchy correctly with expected properties';
$schema->ensure_sequencing_hierarchy(lane => 'esh_lane1', library => 'esh_library1', sample => 'esh_sample2', study => 'esh_study1', group => 'esh_group1', taxon => 'esh_taxon1');
is_deeply [sort { $a <=> $b } map { $_->node_id } $eshlane->related(incoming => { max_depth => 10 })], [$eshlib->node_id, $eshsam->node_id, $eshstu->node_id, $eshgro->node_id, $eshtax->node_id], 'ensure_sequencing_hierarchy() did not change anything when a different sample was supplied';
my $eshsam2 = $schema->get('Sample', { name => 'esh_sample2' });
is $eshsam2, undef, 'and the new sample was not created';
$schema->ensure_sequencing_hierarchy(lane => 'esh_lane1', library => 'esh_library1', sample => 'esh_sample2', study => 'esh_study1', group => 'esh_group1', taxon => 'esh_taxon1', enforce => 1);
$eshsam2 = $schema->get('Sample', { name => 'esh_sample2' });
is_deeply [sort { $a <=> $b } map { $_->node_id } $eshlane->related(incoming => { max_depth => 10 })], [$eshlib->node_id, $eshstu->node_id, $eshgro->node_id, $eshtax->node_id, $eshsam2->node_id], 'ensure_sequencing_hierarchy(enforce => 1) DID change the hierarchy when a different sample was supplied';

# test add() in enqueue mode, and that it maintains history as normal
ok !$schema->add('Sample', { name => 'enqueue1', public_name => 'enqueue1_public_a' }, enqueue => 1), 'add(enqueue => 1) returns nothing';
ok !$schema->get('Sample', { name => 'enqueue1' }), 'it did not add a node to the database';
ok !$schema->add('Sample', { name => 'enqueue2' }, enqueue => 1), 'add(enqueue => 1) worked again';
ok my @queued = $schema->dispatch_queue(), 'could call dispatch_queue()';
is_deeply [sort map { $_->name() } @queued], ['enqueue1', 'enqueue2'], 'dispatch_queue() returned the enqueued nodes';
my @expected_queue = ($schema->get('Sample', { name => 'enqueue1' }), $schema->get('Sample', { name => 'enqueue2' }));
is_deeply [sort map { $_->name() } @expected_queue], ['enqueue1', 'enqueue2'], 'after dispatch_queue() the enqueued nodes are in the database';
$schema->add('Sample', { name => 'enqueue1', public_name => 'enqueue1_public_b' }, enqueue => 1);
$schema->add('Sample', { name => 'enqueue2' }, enqueue => 1);
$schema->add('Sample', { name => 'enqueue3' }, enqueue => 1);
$schema->dispatch_queue();
my @second_queue = ($schema->get('Sample', { name => 'enqueue1' }), $schema->get('Sample', { name => 'enqueue2' }), $schema->get('Sample', { name => 'enqueue3' }));
is_deeply [[sort map { $_->name() } @second_queue], [sort map { $_->node_id() } @second_queue], $second_queue[0]->public_name], [['enqueue1', 'enqueue2', 'enqueue3'], [$expected_queue[0]->node_id, $expected_queue[1]->node_id, $second_queue[2]->node_id], 'enqueue1_public_b'], 'dispatch_queue() was able to update, leave alone and add a new node';
@history = $second_queue[0]->property_history();
is_deeply [$history[0]->{properties}->{public_name}, $history[1]->{properties}->{public_name}], ['enqueue1_public_b', 'enqueue1_public_a'], 'history was maintained on a node updated via dispatch_queue()';

# test search()
my @samples = $schema->search('Sample', { name => '.+queue\d+' });
is scalar(@samples), 3, 'search() worked';
is_deeply [sort map { $_->name() } @samples], ['enqueue1', 'enqueue2', 'enqueue3'], 'search() returned the correct blessed objects';

# before creating the VRPipe schema, test that related() on a VRTrack node can
# return fully-functional VRPipe nodes
$graph->add_schema(namespace => 'VRPipe', label => 'FileSystemElement', unique => [qw(uuid)], indexed => [qw(basename md5 path)]);
my $fse_node = $graph->add_node(namespace => 'VRPipe', label => 'FileSystemElement', properties => { uuid => 'myuuid', basename => 'basename.txt', foo => 'bar' });
$eshsam2->relate_to($fse_node, 'file');
ok my ($rel) = $eshsam2->related(outgoing => { type => 'file' }), 'got the FSE node related to a Sample node before creating a VRPipe schema';
ok $rel->can("basename"), 'the FSE node is fully functional, with a basename method';
is $rel->property('foo'), 'bar', 'and the property() method works as well, useful for allow_anything properties';

# test some VRPipe-specific things
my $vrpipe = VRPipe::Schema->create('VRPipe');
my @paths  = ('/foo/a/cat.txt', '/foo/b/cat.txt', '/foo/a/dog.txt');
my @files  = $vrpipe->get_or_store_filesystem_paths(\@paths);
my %got    = map { $vrpipe->filesystemelement_to_path($_) => 1 } @files;
is_deeply \%got, { '/foo/a/cat.txt' => 1, '/foo/b/cat.txt' => 1, '/foo/a/dog.txt' => 1 }, 'get_or_store_file_paths() and file_to_path() worked correctly';
my $dog_uuid = $files[2]->uuid;
my $llama = $vrpipe->move_filesystemelement($files[2], '/foo/b/llama.txt');
is_deeply [$llama->uuid, $vrpipe->filesystemelement_to_path($llama)], [$dog_uuid, '/foo/b/llama.txt'], 'move_filesystemelement() worked';
$vrpipe->move_filesystemelement('/foo/b', '/foo/a/b');
my $bcat = $vrpipe->get('FileSystemElement', { uuid => $files[1]->uuid });
is $vrpipe->filesystemelement_to_path($bcat), '/foo/a/b/cat.txt', 'move_filesystemelement() can be used to move directories, and it also moves all contained files';
ok my $file = $vrpipe->add('File', { path => '/bar/horse.txt' }), 'adding File is possible';
isa_ok($file, 'VRPipe::Schema::VRPipe::FileSystemElement');
is $file->basename, 'horse.txt', 'File, which is actually a FileSystemElement, has the correct basename';
ok my $gotten_file = $vrpipe->get('File', { path => '/bar/horse.txt' }), 'getting a File is possible';
is $gotten_file->node_id, $file->node_id, 'gotten and added filesystemelement nodes match';
is $vrpipe->get('File', { path => '/bar/horse' }), undef, 'getting a non-existing file does not create it in db';
is $gotten_file->path, '/bar/horse.txt', 'there is a working path method on FileSystemElements';
ok my $new_location = $gotten_file->move('/zar/horse.txt'), 'there is also a working move method';
is_deeply [$new_location->node_id, $new_location->path], [$gotten_file->node_id, '/zar/horse.txt'], 'the move worked correctly';
my $tmp_dir     = $vrpipe->tempdir();
my $source_path = file($tmp_dir, 'source')->stringify;
my $vrfile      = VRPipe::File->create(path => $source_path);
my $fh          = $vrfile->openw;
print $fh "source\n";
$vrfile->close;
my $sym_path   = file($tmp_dir, 'sym')->stringify;
my $cp_path    = file($tmp_dir, 'copy')->stringify;
my $symcp_path = file($tmp_dir, 'symcopy')->stringify;
my $graph_file = $vrpipe->path_to_filesystemelement($source_path);
my $vrsym = VRPipe::File->create(path => $sym_path);
ok $vrfile->symlink($vrsym), 'created a symlink of a file';
ok $vrfile->copy(VRPipe::File->create(path => $cp_path)), 'created a copy of a file';
ok $vrsym->copy(VRPipe::File->create(path => $symcp_path)), 'created a copy of a symlink';
@related = $graph_file->related(outgoing => { max_depth => 2, namespace => 'VRPipe', label => 'FileSystemElement', type => 'symlink|copy' });
my %expected = map { $_->path => 1 } @related;
is_deeply \%expected, { $sym_path => 1, $cp_path => 1, $symcp_path => 1 }, 'there are corresponding nodes in the graph related to the source filesystemelement node';
is $vrpipe->parent_filesystemelement($sym_path)->node_id,   $graph_file->node_id, 'parent_filesystemelement() worked on a symlink path';
is $vrpipe->parent_filesystemelement($cp_path)->node_id,    $graph_file->node_id, 'parent_filesystemelement() worked on a copy path';
is $vrpipe->parent_filesystemelement($symcp_path)->node_id, $graph_file->node_id, 'parent_filesystemelement() worked on a symlink copy path';
my $tmp_sub_dir = dir($tmp_dir, 'sub');
$vrpipe->make_path($tmp_sub_dir);
my $mv_path = file($tmp_sub_dir, 'move')->stringify;
ok $vrfile->move(VRPipe::File->create(path => $mv_path)), 'moved a file';
is_deeply [$vrpipe->parent_filesystemelement($sym_path)->node_id, $vrpipe->parent_filesystemelement($cp_path)->node_id, $vrpipe->parent_filesystemelement($symcp_path)->node_id, $vrpipe->parent_filesystemelement($sym_path)->path, $graph_file->path], [$graph_file->node_id, $graph_file->node_id, $graph_file->node_id, $mv_path, $mv_path], 'after moving the source file the graph of source and symlink and copy is still correct';
my $lrpath     = '/a/path/that/could/be/local/or/remote.txt';
my $local_file = $vrpipe->add('File', { path => $lrpath });
my $irods_file = $vrpipe->add('File', { path => $lrpath, protocol => 'irods:' });
ok $gotten_file = $vrpipe->get('File', { path => $irods_file->protocolless_path, protocol => $irods_file->protocol }), 'you can get() an irods file using its protocol and protocolless_path';
isnt $local_file->uuid, $irods_file->uuid, 'local and remote files with the same absolute paths have different uuids';
my $lf = $vrpipe->get('File', { path => $lrpath });
is $lf->uuid, $local_file->uuid, 'you can get a local file that shares the same path with a remote file';
my $if = $vrpipe->get('File', { path => $lrpath, protocol => 'irods:' });
is $if->uuid, $irods_file->uuid, 'you can get a remote file that shares the same path with a local file';
is $local_file->path, $lrpath, 'path of an local file according to path() looks normal';
is $irods_file->path, "irods:$lrpath", 'path of an irods file according to path() contains the protocol';
is $irods_file->{properties}->{path}, $lrpath, 'path property of the irods file does not contain the protocol';
$vrpipe->move_filesystemelement($irods_file, '/moved/irods/file.txt', protocol => 'irods:');
is $irods_file->path, "irods:/moved/irods/file.txt", 'you can move an irods file';
my $ftp_file = $vrpipe->add('File', { path => $lrpath, protocol => 'ftp://user:pass@host:port' });
is $ftp_file->path, 'ftp://user:pass@host:port' . $lrpath, 'you can specify protocols with passwords and get it back via path()';
is $ftp_file->path(1), $lrpath, 'path(1) returns the path without the protocol';
is $ftp_file->protocolless_path, $lrpath, 'as does protocolless_path()';
my ($ftp_root, $ftp_first_dir) = reverse($ftp_file->related(incoming => { max_depth => 500, namespace => 'VRPipe', label => 'FileSystemElement', type => 'contains' }));
isnt $ftp_root->basename(), 'ftp://user:pass@host:port/', 'however the password is not stored as plain text in the graph db';
my $ftp_file_uuid = $ftp_file->uuid;
$vrpipe->move_filesystemelement($ftp_first_dir, $ftp_first_dir->path(1), protocol => 'ftp://user:changedpass@host:port');
my $ff = $vrpipe->get('File', { path => $lrpath, protocol => 'ftp://user:changedpass@host:port' });
is $ff->uuid, $ftp_file_uuid, 'if your ftp password changed, you could just move the first directory and everything would update automatically';
is $ff->protocol, 'ftp://user:changedpass@host:port', 'protocol() works for an ftp file';
is $ff->protocol(1), 'ftp', 'protocol(1) works for an ftp file';
is $local_file->protocol, 'file:/', 'protocol() works for a local file';
ok $gotten_file = $vrpipe->get('File', { path => $local_file->protocolless_path, protocol => $local_file->protocol }), 'you can get() a local file using its protocol and protocolless_path';
is $local_file->protocol(1), 'file', 'protocol(1) works for a local file';
my $real_local_path = file(qw(t data file.txt))->absolute->stringify;
my $real_local_file = $vrpipe->add('File', { path => $real_local_path });
is $real_local_file->cat_cmd, "cat $real_local_path", 'cat_cmd() works on a local file';
my @expected_data_file_lines = ("a text file\n", "with two lines\n");
ok $fh = $real_local_file->openr, 'openr() method worked on a local file';
my @lines = <$fh>;
is_deeply \@lines, \@expected_data_file_lines, 'the filehandle really works on a local file and let us read the file';
ok $real_local_file->close, 'close() worked on a local file';
SKIP: {
    my $num_tests = 4;
    skip "author-only tests for reading a file from irods", $num_tests unless ($ENV{VRPIPE_AUTHOR_TESTS} && $ENV{VRPIPE_IRODS_TEST_ROOT} && $ENV{VRPIPE_IRODS_TEST_RESOURCE});
    my $irods_root     = $ENV{VRPIPE_IRODS_TEST_ROOT};
    my $irods_resource = $ENV{VRPIPE_IRODS_TEST_RESOURCE};
    system("irm -fr $irods_root > /dev/null 2> /dev/null");
    system("imkdir -p $irods_root");
    system("iput -R $irods_resource $real_local_path $irods_root");
    
    my $real_irods_path = "$irods_root/file.txt";
    my $real_irods_file = $vrpipe->add('File', { path => $real_irods_path, protocol => 'irods:' });
    is $real_irods_file->cat_cmd, "iget $real_irods_path -", 'cat_cmd() works on an irods file';
    ok $fh = $real_irods_file->openr, 'openr() method worked on an irods file';
    @lines = <$fh>;
    is_deeply \@lines, \@expected_data_file_lines, 'the filehandle really works on an irods file and let us read the file';
    ok $real_irods_file->close, 'close() worked on an irods file';
    
    system("irm -fr $irods_root");
}

# make a traditional mysql vrpipe StepState so we can test that
# ensure_state_hierarchy() can represent the same thing in the graph database
my $pipeline = VRPipe::Pipeline->create(name => 'test_pipeline');
my $datasource = VRPipe::DataSource->create(type => 'fofn_with_metadata', method => 'grouped_by_metadata', options => { metadata_keys => 'study|sample' }, source => file(qw(t data datasource.fofn_with_metadata))->absolute);
my $setup = VRPipe::PipelineSetup->create(name => 'ps1', datasource => $datasource, output_root => '/tmp/out', pipeline => $pipeline, active => 0);
$datasource->elements;
my $ss1 = VRPipe::StepState->create(dataelement => 1, pipelinesetup => 1, stepmember => 1);
my $ss2 = VRPipe::StepState->create(dataelement => 2, pipelinesetup => 1, stepmember => 1);
ok my $graph_ss = $vrpipe->ensure_state_hierarchy($ss1), 'ensure_state_hierarchy() worked';
$graph_ss = $vrpipe->ensure_state_hierarchy($ss2);
ok $graph_ss = $vrpipe->ensure_state_hierarchy($ss2), 'ensure_state_hierarchy() worked twice on the same StepState';
ok my $ps = $vrpipe->get('PipelineSetup', { name => 'ps1' }), 'a PipelineSetup was created in the graph database';
@related = $ps->related(outgoing => { max_depth => 500 });
is scalar(@related), 21, 'all the related nodes were also created';
#*** more detailed, specific tests to make sure the graph is correct? checked visually it's fine...

# test VRTrack schema's add_file method, which passes through to VRPipe schema
my $bar_snake = '/bar/snake.txt';
ok my $vrtrack_file = $schema->add_file($bar_snake), 'add_file() method on the VRTrack schema worked';
is $vrtrack_file->path, $bar_snake, 'the path() method worked on what that returned';
ok my $vrtrack_irods_file = $schema->add_file($bar_snake, 'irods:'), 'add_file() method on the VRTrack schema worked with a protocol';
isnt $vrtrack_file->uuid, $vrtrack_irods_file->uuid, 'local and irods files with the same path are different nodes';
is $schema->get_file($bar_snake)->uuid, $vrtrack_file->uuid, 'VRTrack get_file gets the correct file with no protocol';
is $schema->get_file($bar_snake, 'irods:')->uuid, $vrtrack_irods_file->uuid, 'VRTrack get_file gets the correct file with a protocol';
is $vrtrack_irods_file->path, 'irods:' . $bar_snake, 'path() is correct for an irods file';

# check there are no issues supplying ints versus strings to unique values
$schema->add('Study', { id => 2625,   name => 'str_vs_int_test', accession => 'svitacc' });
$schema->add('Study', { id => '2625', name => 'str_vs_int_test', accession => 'svitacc' });
my @studies = $schema->get("Study", { name => 'str_vs_int_test' });
is scalar(@studies), 1, 'providing int or string as the unique when adding a node only results in 1 new node';
my $s_via_int = $schema->get("Study", { id => 2625 });
my $s_via_str = $schema->get("Study", { id => '2625' });
is_deeply [$s_via_int ? $s_via_int->{id} : 0, $s_via_str ? $s_via_str->{id} : 0], [$studies[0]->{id}, $studies[0]->{id}], 'a study can be gotten via its unique id as an int or a string';
# also check there are no issues supplying strings with colons in them
$schema->add('Study', { id => 2626, name => 'Study ERP006001: Deep sequencing of HGDP samples on the Illumina X10', accession => 'colon' });
$schema->add('Study', { id => 2626, name => 'Study ERP006001: Deep sequencing of HGDP samples on the Illumina X10', accession => 'colon' });
my $study = $schema->get("Study", { id => 2626 });
is_deeply $study->{properties}, { id => 2626, name => 'Study ERP006001: Deep sequencing of HGDP samples on the Illumina X10', accession => 'colon' }, 'colons in strings do not break nodes';

# check that fixes for the above didn't break cypher queries with node ids
my $donor = $schema->add('Donor', { id => 'd1' }, outgoing => { type => 'has', node => $sample });
ok my $extra_info_node = $schema->get_node_by_id_with_extra_info('Donor', $donor->{id}), 'cypher queries with node ids work';

# test divorce_from()
my $samson  = $schema->add('Group', { name => 'Samson' });
my $delilah = $schema->add('Group', { name => 'Delilah' });
$samson->relate_to($delilah, 'loves');
$delilah->relate_to($samson, 'betrays');
@related = $samson->related(outgoing => {});
is_deeply [sort map { $_->node_id } @related], [$delilah->node_id], 'Samson is related to Delilah';
@related = $delilah->related(outgoing => {});
is_deeply [sort map { $_->node_id } @related], [$samson->node_id], 'Delilah is related to Samson';
$delilah->divorce_from($samson, 'loves');
@related = $samson->related(outgoing => {});
is scalar(@related), 0, 'Samson is no longer related to Delilah after divorce_from on loves';
@related = $delilah->related(outgoing => {});
is_deeply [sort map { $_->node_id } @related], [$samson->node_id], 'Delilah is still related to Samson';
$samson->divorce_from($delilah);
@related = $samson->related(outgoing => {});
is scalar(@related), 0, 'Samson is sill not related to Delilah after divorce_from on all';
@related = $delilah->related(outgoing => {});
is scalar(@related), 0, 'and now Delilah is not related to Samson';

# check that we can have arrayref properties
ok my $ganode = $schema->add('Group', { name => 'foo', qc_fail_reasons => ['one', 'two', 'three'] }), 'we can set a property with an array while creating the node';
is_deeply $ganode->qc_fail_reasons(), ['one', 'two', 'three'], 'we can get array values of a property';
ok $ganode->qc_fail_reasons(['four', 'five', 'six']), 'we can set a property with an array using the auto get/setter';
$ganode = $schema->get('Group', { name => 'foo' });
is_deeply $ganode->qc_fail_reasons(), ['four', 'five', 'six'], 'setting an array with the get/setter really worked';

# test removing properties, and how history works with that
my $bps = $schema->add('Sample', { name => 'bps', public_name => ['b', 'p', 's'] });
throws_ok { $bps->remove_property('name') } qr/Property 'name' supplied, but that's unique for schema VRTrack::Sample and can't be changed/, 'remove_property() throws when given a unique property';
ok $bps->remove_property('public_name'), 'remove_property() worked on a normal property';
is_deeply [$bps->{properties}, $bps->changed()], [{ name => 'bps' }, { public_name => [['b', 'p', 's'], undef] }], 'properties and changed() after removal of a property give expected results';
$bps->add_properties({ public_name => 'pn' });
@history = $bps->property_history('public_name');
@history_properties = map { $_->{properties} } @history;
is_deeply [@history_properties], [{ public_name => 'pn' }, {}, { public_name => ['b', 'p', 's'] }], 'property_history() works correctly on a property that had at one point been removed';

# test get_sequencing_hierarchy; first create more hierarchies, some of which
# have samples belonging to more than 1 study
$schema->ensure_sequencing_hierarchy(lane => 'esh_lane2', library => 'esh_library2', sample => 'esh_sample2', study => 'esh_study1', group => 'esh_group1', taxon => 'esh_taxon1');
$hierarchy = $schema->ensure_sequencing_hierarchy(lane => 'esh_lane3', library => 'esh_library3', sample => 'esh_sample2', study => 'esh_study2', group => 'esh_group1', taxon => 'esh_taxon1');
$schema->ensure_sequencing_hierarchy(lane => 'esh_lane4', library => 'esh_library4', sample => 'esh_sample2', study => 'esh_study3', group => 'esh_group1', taxon => 'esh_taxon1');
($rel) = @{ $graph->related($hierarchy->{study}, undef, undef, {})->{relationships} };
$graph->relationship_set_properties($rel, { preferred => 1 });
my $lane3_file = $vrpipe->add('File', { path => '/seq/3/3.bam' });
$hierarchy->{lane}->relate_to($lane3_file, 'aligned');
ok $hierarchy = $schema->get_sequencing_hierarchy($lane3_file), 'get_sequencing_hierarchy() seemed to work';
%hierarchy_props = map { $_ => $hierarchy->{$_}->properties } keys %{$hierarchy};
is_deeply \%hierarchy_props, { lane => { unique => 'esh_lane3', lane => 'esh_lane3' }, library => { id => 'esh_library3' }, sample => { name => 'esh_sample2' }, study => { id => 'esh_study2' }, taxon => { id => 'esh_taxon1' } }, 'get_sequencing_hierarchy() returned the correct nodes for a bam file, including the preferred study';
$schema->ensure_sequencing_hierarchy(lane => 'esh_lane5', library => 'esh_library5', sample => 'esh_sample5', study => 'esh_study5', group => 'esh_group2', taxon => 'esh_taxon1');
$hierarchy = $schema->ensure_sequencing_hierarchy(lane => 'esh_lane6', library => 'esh_library6', sample => 'esh_sample5', study => 'esh_study6', group => 'esh_group2', taxon => 'esh_taxon1');
$schema->ensure_sequencing_hierarchy(lane => 'esh_lane7', library => 'esh_library7', sample => 'esh_sample5', study => 'esh_study7', group => 'esh_group2', taxon => 'esh_taxon1');
my $lane6_file = $vrpipe->add('File', { path => '/seq/6/6.bam' });
$hierarchy->{lane}->relate_to($lane6_file, 'aligned');
$hierarchy = $schema->get_sequencing_hierarchy($lane6_file);
my $hstudy = delete $hierarchy->{study};
%hierarchy_props = map { $_ => $hierarchy->{$_}->properties } keys %{$hierarchy};
is_deeply \%hierarchy_props, { lane => { unique => 'esh_lane6', lane => 'esh_lane6' }, library => { id => 'esh_library6' }, sample => { name => 'esh_sample5' }, taxon => { id => 'esh_taxon1' } }, 'get_sequencing_hierarchy() returned the correct nodes for another bam file';
like $hstudy->id, qr/esh_study[567]/, 'since no preferred study was specified, a random one was returned';
$hierarchy = $schema->ensure_sequencing_hierarchy(lane => 'esh_lane8', library => 'esh_library8', sample => 'esh_sample8', study => 'esh_study8', group => 'esh_group3', taxon => 'esh_taxon1');
my $lane8_file = $vrpipe->add('File', { path => '/seq/8/8.bam' });
$hierarchy->{lane}->relate_to($lane8_file, 'aligned');
$hierarchy = $schema->get_sequencing_hierarchy($lane8_file);
%hierarchy_props = map { $_ => $hierarchy->{$_}->properties } keys %{$hierarchy};
is_deeply \%hierarchy_props, { lane => { unique => 'esh_lane8', lane => 'esh_lane8' }, library => { id => 'esh_library8' }, sample => { name => 'esh_sample8' }, study => { id => 'esh_study8' }, taxon => { id => 'esh_taxon1' } }, 'get_sequencing_hierarchy() returned the correct nodes for another bam file, in the simple case of there being only 1 study for the sample';

my $hierarchy_meta = $schema->node_and_hierarchy_properties($lane8_file);
is_deeply $hierarchy_meta, { vrtrack_lane_unique => 'esh_lane8', vrtrack_lane_lane => 'esh_lane8', vrtrack_library_id => 'esh_library8', vrtrack_sample_name => 'esh_sample8', vrtrack_study_id => 'esh_study8', vrtrack_taxon_id => 'esh_taxon1' }, 'node_and_hierarchy_properties() worked';

# test file_qc_nodes() by first manually adding the nodes that
# npg_cram_stats_parser step adds
my $qc_file = $vrpipe->add('File', { path => '/seq/8/8.stats' });
$lane8_file->relate_to($qc_file, 'qc_file');
$schema->add(
    'Bam_Stats',
    {
        mode                  => 'normal',
        options               => '-opt',
        date                  => 1443015733,
        'reads QC failed'     => 99,
        'raw total sequences' => 10001
    },
    incoming => { type => 'summary_stats', node => $qc_file }
);
$qc_file = $vrpipe->add('File', { path => '/seq/8/8.genotype.json' });
$lane8_file->relate_to($qc_file, 'qc_file');
$schema->add(
    'Genotype',
    {
        date                 => 1443015733,
        pass                 => 1,
        expected_sample_name => 'foo',
        matched_sample_name  => 'foo'
    },
    incoming => { type => 'genotype_data', node => $qc_file }
);
$qc_file = $vrpipe->add('File', { path => '/seq/qc/8/8.verify_bam_id.json' });
$lane8_file->relate_to($qc_file, 'qc_file');
$schema->add(
    'Verify_Bam_ID',
    {
        date    => 1443015733,
        pass    => 0,
        freemix => 'foo'
    },
    incoming => { type => 'verify_bam_id_data', node => $qc_file }
);

my $qc_nodes = $schema->file_qc_nodes($lane8_file);
my %qc_props = map { $_ => $qc_nodes->{$_}->properties } keys %{$qc_nodes};
delete $qc_props{bam_stats}->{uuid};
delete $qc_props{genotype}->{uuid};
delete $qc_props{verify_bam_id}->{uuid};
is_deeply \%qc_props, { bam_stats => { mode => 'normal', options => '-opt', date => 1443015733, 'reads QC failed' => 99, 'raw total sequences' => 10001 }, genotype => { date => 1443015733, pass => 1, expected_sample_name => 'foo', matched_sample_name => 'foo' }, verify_bam_id => { date => 1443015733, pass => 0, freemix => 'foo' } }, 'file_qc_nodes() returned the correct nodes';

my $vrtrack_meta = $schema->vrtrack_metadata($lane8_file);
is_deeply $vrtrack_meta, { %{$hierarchy_meta}, vrtrack_bam_stats_mode => 'normal', vrtrack_bam_stats_options => '-opt', vrtrack_bam_stats_date => 1443015733, 'vrtrack_bam_stats_reads QC failed' => 99, 'vrtrack_bam_stats_raw total sequences' => 10001, vrtrack_genotype_date => 1443015733, vrtrack_genotype_pass => 1, vrtrack_genotype_expected_sample_name => 'foo', vrtrack_genotype_matched_sample_name => 'foo', vrtrack_verify_bam_id_date => 1443015733, vrtrack_verify_bam_id_pass => 0, vrtrack_verify_bam_id_freemix => 'foo' }, 'vrtrack_metadata() worked';

# also test the pass-through from VRPipe::File to get vrtrack metadata
$vrfile = VRPipe::File->create(path => '/seq/8/8.bam');
$vrfile->add_metadata({ sql_meta => 'sql_value' });
is_deeply $vrfile->metadata(undef, include_vrtrack => 1), { %{$vrtrack_meta}, sql_meta => 'sql_value' }, 'VRPipe::File->metadata(undef, include_vrtrack => 1) worked';

exit;
