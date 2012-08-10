#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 49;
    use VRPipeTest;
    
    use_ok('VRPipe::DataSourceFactory');
}

# list
ok my $ds = VRPipe::DataSource->create(type    => 'list',
                                       method  => 'all',
                                       source  => file(qw(t data datasource.list))->absolute->stringify,
                                       options => {}),
  'could create a list datasource';

my @results;
foreach my $element (@{ get_elements($ds) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ line => 'foo' }, { line => 'bar' }, { line => 'henry' }], 'got all results correctly';

$ds = VRPipe::DataSource->create(type    => 'list',
                                 method  => 'all',
                                 source  => file(qw(t data datasource.list))->absolute->stringify,
                                 options => { skip_comments => 0 });

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ line => 'foo' }, { line => 'bar' }, { line => '# comment' }, { line => 'henry' }], 'got even more results with extra options';

# fofn
ok $ds = VRPipe::DataSource->create(type    => 'fofn',
                                    method  => 'all',
                                    source  => file(qw(t data datasource.fofn))->absolute->stringify,
                                    options => {}),
  'could create a fofn datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => [file('t', 'data', 'file.bam')->absolute] }, { paths => [file('t', 'data', 'file.cat')->absolute] }, { paths => [file('t', 'data', 'file.txt')->absolute] }], 'got correct results for fofn all';

ok $ds = VRPipe::DataSource->create(type    => 'fofn',
                                    method  => 'group_all',
                                    source  => file(qw(t data datasource.fofn))->absolute->stringify,
                                    options => {}),
  'could create a fofn datasource with group_all method';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => [file('t', 'data', 'file.bam')->absolute, file('t', 'data', 'file.cat')->absolute, file('t', 'data', 'file.txt')->absolute] }], 'got correct results for fofn group_all';

# delimited
ok $ds = VRPipe::DataSource->create(type    => 'delimited',
                                    method  => 'grouped_single_column',
                                    source  => file(qw(t data datasource.fastqs))->absolute->stringify,
                                    options => { delimiter => "\t",
                                                 group_by  => 1,
                                                 column    => 2 }),
  'could create a delimited datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => [file(qw(t data 2822_6_1.fastq))->absolute, file(qw(t data 2822_6_2.fastq))->absolute], group => '2822_6' }, { paths => [file(qw(t data 2822_7_1.fastq))->absolute, file(qw(t data 2822_7_2.fastq))->absolute], group => '2822_7' }, { paths => [file(qw(t data 2823_4_1.fastq))->absolute, file(qw(t data 2823_4_2.fastq))->absolute], group => '2823_4' }, { paths => [file(qw(t data 8324_8_1.fastq))->absolute, file(qw(t data 8324_8_2.fastq))->absolute], group => '8324_8' }], 'got correct results for delimited grouped_single_column';

# delimited all columns
ok $ds = VRPipe::DataSource->create(type    => 'delimited',
                                    method  => 'all_columns',
                                    source  => file(qw(t data datasource.2col))->absolute->stringify,
                                    options => { delimiter => "\t" }),
  'could create a delimited datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => [file(qw(t data file.txt))->absolute, file(qw(t data file2.txt))->absolute] }, { paths => [file(qw(t data file3.txt))->absolute, file(qw(t data file.cat))->absolute] }], 'got correct results for delimited all_columns';

# test get_methods and method_options
$ds = $ds->_source_instance;
is_deeply [sort $ds->get_methods], [qw(all all_columns grouped_single_column single_column)], 'get_methods returned all the expected methods';
is_deeply [$ds->method_options('single_column')], [['named', 'delimiter', 1, undef, 'Str'], ['named', 'column', 1, undef, 'PositiveInt'], ['named', 'column_is_path', 0, 1, 'Bool']], 'method_options showed us what options the single_column method takes';

# fofn_with_metadata
my @fwm_paths = ('/a/path/7816_3#95.bam', '/a/path/7413_5#95.bam', '/a/path/8312_5#95.bam');
my %fwm_common_meta = (center_name => 'SC', study => 'ERP000979', platform => 'ILLUMINA');
VRPipe::File->create(path => $fwm_paths[0])->add_metadata({ library => 'foo' });
ok $ds = VRPipe::DataSource->create(type    => 'fofn_with_metadata',
                                    method  => 'all',
                                    source  => file(qw(t data datasource.fofn_with_metadata))->absolute->stringify,
                                    options => {}),
  'could create a fofn_with_metadata datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, $element->result);
}
is_deeply [@results, VRPipe::File->get(path => $fwm_paths[0])->metadata, VRPipe::File->get(path => $fwm_paths[1])->metadata, VRPipe::File->get(path => $fwm_paths[2])->metadata], [{ paths => [$fwm_paths[0]] }, { paths => [$fwm_paths[1]] }, { paths => [$fwm_paths[2]] }, { %fwm_common_meta, sample => 'JB953', library => '4858080', lane => '7816_3#95' }, { %fwm_common_meta, sample => 'JB953', library => '4074406', lane => '7413_5#95' }, { %fwm_common_meta, sample => 'JB951', library => '4074399', lane => '8312_5#95' }], 'got correct results for fofn_with_metadata all, and the metadata on the files was correct';

ok $ds = VRPipe::DataSource->create(type    => 'fofn_with_metadata',
                                    method  => 'group_all',
                                    source  => file(qw(t data datasource.fofn_with_metadata))->absolute->stringify,
                                    options => {}),
  'could create a fofn_with_metadata datasource with group_all method';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, $element->result);
}
is_deeply [@results, VRPipe::File->get(path => $fwm_paths[0])->metadata, VRPipe::File->get(path => $fwm_paths[1])->metadata, VRPipe::File->get(path => $fwm_paths[2])->metadata], [{ paths => [$fwm_paths[0], $fwm_paths[1], $fwm_paths[2]] }, { %fwm_common_meta, sample => 'JB953', library => '4858080', lane => '7816_3#95' }, { %fwm_common_meta, sample => 'JB953', library => '4074406', lane => '7413_5#95' }, { %fwm_common_meta, sample => 'JB951', library => '4074399', lane => '8312_5#95' }], 'got correct results for fofn_with_metadata group_all, and the metadata on the files was correct';

ok $ds = VRPipe::DataSource->create(type    => 'fofn_with_metadata',
                                    method  => 'grouped_by_metadata',
                                    source  => file(qw(t data datasource.fofn_with_metadata))->absolute->stringify,
                                    options => { metadata_keys => 'study|sample' }),
  'could create a fofn_with_metadata grouped_by_metadata datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, $element->result);
}
is_deeply [sort { $a->{group} cmp $b->{group} } @results], [{ paths => [$fwm_paths[2]], group => 'ERP000979|JB951' }, { paths => [$fwm_paths[0], $fwm_paths[1]], group => 'ERP000979|JB953' }], 'got correct results for fofn_with_metadata grouped_by_metadata';

# sequence_index
ok $ds = VRPipe::DataSource->create(type    => 'sequence_index',
                                    method  => 'lane_fastqs',
                                    source  => file(qw(t data datasource.sequence_index))->absolute->stringify,
                                    options => { local_root_dir => dir('./')->absolute->stringify }),
  'could create a sequence_index datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, $element->result);
}
is_deeply \@results, [{ paths => [file(qw(t data 2822_6.fastq))->absolute, file(qw(t data 2822_6_1.fastq))->absolute, file(qw(t data 2822_6_2.fastq))->absolute], lane => '2822_6' }, { paths => [file(qw(t data 2822_7_1.fastq))->absolute, file(qw(t data 2822_7_2.fastq))->absolute], lane => '2822_7' }, { paths => [file(qw(t data 2823_4_1.fastq))->absolute, file(qw(t data 2823_4_2.fastq))->absolute], lane => '2823_4' }, { paths => [file(qw(t data 8324_8_1.fastq))->absolute, file(qw(t data 8324_8_2.fastq))->absolute], lane => '8324_8' }], 'got correct results for sequence_index lane_fastqs';
my $vrfile = VRPipe::File->get(path => file(qw(t data 2822_6_1.fastq))->absolute);
my $meta = $vrfile->metadata;
is_deeply $meta,
  { expected_md5   => 'f1826489facca0d0bdf02d9586b493f6',
    lane           => '2822_6',
    study          => 'STUDY01',
    study_name     => 'my study name',
    center_name    => 'SC',
    sample_id      => 'SAMPLEID01',
    sample         => 'SAMPLE01',
    population     => 'POP',
    platform       => 'ILLUMINA',
    library        => 'LIB01',
    insert_size    => 200,
    withdrawn      => 0,
    reads          => 200,
    bases          => 12200,
    analysis_group => 'low coverage',
    paired         => 1,
    mate           => file(qw(t data 2822_6_2.fastq))->absolute->stringify },
  'a VRPipe::File created by source has the correct metadata';

# test a special vrtrack test database; these tests are meant for the author
# only, but will also work for anyone with a working VRTrack::Factory setup
SKIP: {
    my $num_tests = 26;
    skip "author-only tests for a VRTrack datasource", $num_tests unless $ENV{VRPIPE_VRTRACK_TESTDB};
    eval "require VRTrack::Factory;";
    skip "VRTrack::Factory not loading", $num_tests if $@;
    
    # create the vrtrack db
    my %cd  = VRTrack::Factory->connection_details('rw');
    my @sql = VRTrack::VRTrack->schema();
    open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
    print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    foreach my $sql (@sql) {
        print $mysqlfh $sql;
    }
    close($mysqlfh);
    
    # populate it *** uses update_vrmeta.pl, which may not be supplied with
    #                 VRTrack in the future...
    system("update_vrmeta.pl --samples t/data/vrtrack.samples --index t/data/vrtrack.sequence.index --database $ENV{VRPIPE_VRTRACK_TESTDB} > /dev/null 2> /dev/null");
    
    # alter processed on the lanes to enable useful tests
    my $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
    my $lanes = $vrtrack->processed_lane_hnames();
    my %expectations;
    for my $i (1 .. 60) {
        my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lanes->[$i - 1]);
        my $name = $lane->hierarchy_name;
        
        if ($i <= 50) {
            $lane->is_processed('import' => 1);
            push(@{ $expectations{import} }, $name);
            if ($i > 10) {
                $lane->is_processed('qc' => 1);
                push(@{ $expectations{qc} }, $name);
                if ($i > 45) {
                    $lane->qc_status('passed');
                    push(@{ $expectations{qc_status_passed} }, $name);
                }
            }
            if ($i > 20) {
                $lane->is_processed('mapped' => 1);
                push(@{ $expectations{mapped} }, $name);
            }
            if ($i > 30) {
                $lane->is_processed('improved' => 1);
                push(@{ $expectations{improved} }, $name);
            }
            if ($i > 40) {
                $lane->is_processed('stored' => 1);
                push(@{ $expectations{stored} }, $name);
            }
        }
        elsif ($i > 55) {
            $lane->is_withdrawn(1);
            push(@{ $expectations{withdrawn} }, $name);
        }
        
        $lane->update;
    }
    
    ok $ds = VRPipe::DataSource->create(type    => 'vrtrack',
                                        method  => 'lanes',
                                        source  => $ENV{VRPIPE_VRTRACK_TESTDB},
                                        options => { import => 1, mapped => 0 }),
      'could create a vrtrack datasource';
    my $results = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $results++;
    }
    is $results, 20, 'got correct number of results for vrtrack lanes mapped => 0';

    ## tests for _has_changed
    ok(!$ds->_source_instance->_has_changed, 'vrtrack datasource _has_changed gives no change');
    
    # create a new lane
    $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
    my $new_lane = VRTrack::Lane->create($vrtrack, 'new_lane');
    $new_lane->is_withdrawn(0);
    $new_lane->library_id(16);
    $new_lane->update;
    ok($ds->_source_instance->_has_changed, 'vrtrack datasource _has_changed gives change after new lane insertion in test vrtrack db');
    
    # Go back to unchanged state by deleting this lane. Check we don't have any changes
    $new_lane->delete;
    ok(!$ds->_source_instance->_has_changed, 'vrtrack datasource _has_change gives no change after inserted lane deleted in test vrtrack db');
    
    # add a file to another lane and check for changes
    my $lane_to_add_file_for = VRTrack::Lane->new_by_name($vrtrack, 'ERR003040');
    my $newfile = $lane_to_add_file_for->add_file('new.fastq');
    is $ds->_source_instance->_has_changed, 0, 'vrtrack datasource _has_changed gives no change after adding a file in test vrtrack db, with method lanes';
    $newfile->delete;
    
    # change some md5 sums in the files
    my $file = VRTrack::File->new_by_hierarchy_name($vrtrack, 'ERR003038.filt.fastq.gz');
    # check for changes
    $file->md5('34c009157187c5d9a7e976563ec1bad8');
    $file->update;
    is $ds->_source_instance->_has_changed, 0, 'datasource _has_changed got no change after md5 change in file table in test vrtrack db, with method lanes';
    $file->md5('cac33e4fc8ff2801978cfd5a223f5064');
    $file->update;
    
    # if we change something that doesn't indicate a real change as far as our
    # options are concerned, it shouldn't come up as _has_changed
    $lane_to_add_file_for->qc_status('pending');
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 0, 'changing qc_status when it was not an option is ignored';
    $lane_to_add_file_for->is_processed('mapped', 1);
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 1, 'changing is_processed when it was an option is detected';
    $lane_to_add_file_for->is_processed('mapped', 0);
    $lane_to_add_file_for->qc_status('no_qc');
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 0, 'reverting lane back gives no change';
    $lane_to_add_file_for->is_withdrawn(1);
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 1, 'changing withdrawn is always detected';
    $lane_to_add_file_for->is_withdrawn(0);
    $lane_to_add_file_for->update;
    
    $ds = VRPipe::DataSource->create(type    => 'vrtrack',
                                     method  => 'lanes',
                                     source  => $ENV{VRPIPE_VRTRACK_TESTDB},
                                     options => { qc_status => 'pending' });
    $results = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $results++;
    }
    is $results, 0, 'initially got no results for vrtrack lanes qc_status => pending';
    $lane_to_add_file_for->qc_status('pending');
    $lane_to_add_file_for->update;
    $results = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $results++;
    }
    is $results, 1, 'after changing a lane to qc pending, got 1 dataelement';
    
    $ds = VRPipe::DataSource->create(type    => 'vrtrack',
                                     method  => 'lanes',
                                     source  => $ENV{VRPIPE_VRTRACK_TESTDB},
                                     options => {});
    $results = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $results++;
    }
    is $results, 55, 'with no options we get all the active lanes in the db';
    $lane_to_add_file_for->qc_status('pending');
    $lane_to_add_file_for->is_processed('mapped', 1);
    $lane_to_add_file_for->raw_bases(50);
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 0, 'we can change lots of things without being detected';
    $lane_to_add_file_for->is_withdrawn(1);
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 1, 'but withdrawn is still detected';
    
    $lane_to_add_file_for->qc_status('no_qc');
    $lane_to_add_file_for->is_processed('mapped', 0);
    $lane_to_add_file_for->raw_bases(1473337566);
    $lane_to_add_file_for->is_withdrawn(0);
    $lane_to_add_file_for->update;
    
    # lane_fastqs tests
    ok $ds = VRPipe::DataSource->create(type    => 'vrtrack',
                                        method  => 'lane_fastqs',
                                        source  => $ENV{VRPIPE_VRTRACK_TESTDB},
                                        options => { import => 1, mapped => 0, local_root_dir => dir('t')->absolute->stringify, library_regex => 'g1k-sc-NA19190-YRI-1\|SC\|SRP000542\|NA19190' }),
      'could create a vrtrack datasource';
    
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, $element->result);
    }
    
    is_deeply $results[0],
      { paths => [file(qw(t data NA19190 sequence_read ERR003199.filt.fastq.gz))->absolute, file(qw(t data NA19190 sequence_read ERR003199_1.filt.fastq.gz))->absolute, file(qw(t data NA19190 sequence_read ERR003199_2.filt.fastq.gz))->absolute],
        lane  => 'ERR003199' },
      'got correct results for vrtrack lane_fastqs';
    my $vrfile = VRPipe::File->get(path => file(qw(t data NA19190 sequence_read ERR003199.filt.fastq.gz))->absolute);
    my $meta = $vrfile->metadata;
    is_deeply $meta,
      { 'bases'        => '1696783',
        'withdrawn'    => 0,
        'project'      => 'SRP000542',
        'species'      => 'Homo sapiens',
        'population'   => 'YRI',
        'paired'       => 1,
        'reads'        => '45859',
        'library'      => 'g1k_sc_NA19190_YRI_1',
        'lane_id'      => '101',
        'individual'   => 'NA19190',
        'center_name'  => 'SC',
        'sample'       => 'NA19190',
        'platform'     => 'SLX',
        'expected_md5' => 'dfa4364855815d7433c451a87f0520d0',
        'study'        => 'SRP000542',
        'lane'         => 'ERR003199',
        'insert_size'  => 175 },
      'a VRPipe::File created by vrtrack datasource has the correct metadata';
    
    $newfile = $lane_to_add_file_for->add_file('new.fastq');
    $newfile->type(0);
    $newfile->update;
    is $ds->_source_instance->_has_changed, 1, 'vrtrack datasource _has_changed gives change after adding a fastq file in test vrtrack db, with method lane_fastqs';
    $newfile->delete;
    
    $newfile = $lane_to_add_file_for->add_file('new.bam');
    $newfile->type(5);
    $newfile->update;
    is $ds->_source_instance->_has_changed, 0, 'vrtrack datasource _has_changed gives no change after adding a bam file in test vrtrack db, with method lane_fastqs';
    $newfile->delete;
    
    # change some md5 sums in the files
    $file = VRTrack::File->new_by_hierarchy_name($vrtrack, 'ERR003038.filt.fastq.gz');
    # check for changes
    $file->md5('34c009157187c5d9a7e976563ec1bad8');
    $file->update;
    is $ds->_source_instance->_has_changed, 1, 'datasource _has_changed got change after md5 change in file table in test vrtrack db, with method lane_fastqs';
    $file->md5('cac33e4fc8ff2801978cfd5a223f5064');
    $file->update;
    is $ds->_source_instance->_has_changed, 0, 'reverting file md5 back gives no change';
    
    # if we change the insert_size in vrtrack, this should cause the datasource
    # to reset the elements. Critically, the vrfile metadata should get updated,
    # or else we'd be stuck in an infinite loop of always detecting the
    # insert_size change and reseting
    $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
    my $lib = VRTrack::Library->new_by_name($vrtrack, 'g1k-sc-NA19190-YRI-1|SC|SRP000542|NA19190');
    $lib->fragment_size_from(200);
    $lib->fragment_size_to(200);
    $lib->update;
    # (we also add a lane to force datasource to notice a change)
    VRTrack::Lane->create($vrtrack, 'new_lane2');
    get_elements($ds);
    $vrfile->reselect_values_from_db;
    is $vrfile->metadata->{insert_size}, '200', 'changing insert_size in vrtrack changes insert_size metadata on vrpipe files';
    
    # test getting improved bams that passed qc
    $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
    my %expected_qc_passed_improved_bams;
    my %passed_hnames = map { $_ => 1 } @{ $expectations{qc_status_passed} };
    foreach my $hname (@{ $expectations{qc} }) {
        my $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $hname);
        my $improved_bam = VRPipe::File->create(path => file($hname . '.improved.bam')->absolute);
        my $vrfile_name  = 'VRPipe::File::' . $improved_bam->id;
        my $md5          = 'an_md5_' . $improved_bam->id;
        $vrfile = $vrlane->add_file($vrfile_name);
        $vrfile->type(5);
        $vrfile->md5($md5);
        $vrfile->update;
        $improved_bam->md5($md5);
        $improved_bam->update;
        
        if (exists $passed_hnames{$hname}) {
            $expected_qc_passed_improved_bams{ $improved_bam->path->stringify } = 1;
        }
        
        # to give us more than 1 group to test, change the library of one of the
        # lanes
        if ($hname eq 'ERR008838') {
            $vrlane->library_id(1);
            $vrlane->update;
        }
    }
    my %expected_groups = ('SRP000546|SRP000546|NA18633|HUMsgR3AIDCAASE'      => 1,
                           'SRP000547|SRP000547|NA07056|g1k_sc_NA07056_CEU_1' => 1);
    $ds = VRPipe::DataSource->create(type    => 'vrtrack',
                                     method  => 'lane_improved_bams',
                                     source  => $ENV{VRPIPE_VRTRACK_TESTDB},
                                     options => { qc                => 1,
                                                  qc_status         => 'passed',
                                                  group_by_metadata => 'project|study|sample|library' });
    my %actual_qc_passed_improved_bams;
    my %actual_groups;
    foreach my $element (@{ get_elements($ds) }) {
        my $result = $element->result;
        $actual_groups{ $result->{group} || 'no_group' } = 1;
        foreach my $path (@{ $result->{paths} }) {
            $actual_qc_passed_improved_bams{$path} = 1;
        }
    }
    is_deeply \%actual_qc_passed_improved_bams, \%expected_qc_passed_improved_bams, 'an improved bam qc passed datasource gave the expected files';
    is_deeply \%actual_groups, \%expected_groups, 'a group_by_metadata datasource gave the expected groups';
}
exit;
