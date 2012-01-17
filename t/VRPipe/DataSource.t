#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Cwd;

BEGIN {
    use Test::Most tests => 24;
    use VRPipeTest;
    
    use_ok('VRPipe::DataSourceFactory');
}

# list
ok my $ds = VRPipe::DataSource->get(type => 'list',
                                    method => 'all',
                                    source => file(qw(t data datasource.list))->absolute->stringify,
                                    options => {}), 'could create a list datasource';

my @results;
foreach my $element (@{$ds->elements}) {
    push(@results, $element->result);
}
is_deeply \@results, [{line => 'foo'}, {line => 'bar'}, {line => 'henry'}], 'got all results correctly';

$ds = VRPipe::DataSource->get(type => 'list',
                              method => 'all',
                              source => file(qw(t data datasource.list))->absolute->stringify,
                              options => {skip_comments => 0});

@results = ();
foreach my $element (@{$ds->elements}) {
    push(@results, $element->result);
}
is_deeply \@results, [{line => 'foo'}, {line => 'bar'}, {line => '# comment'}, {line => 'henry'}], 'got even more results with extra options';

# fofn
ok $ds = VRPipe::DataSource->get(type => 'fofn',
                                 method => 'all',
                                 source => file(qw(t data datasource.fofn))->absolute->stringify,
                                 options => {}), 'could create a fofn datasource';

@results = ();
foreach my $element (@{$ds->elements}) {
    push(@results, $element->result);
}
my $cwd = cwd();
is_deeply \@results, [{paths => [file($cwd, 't', 'data', 'file.bam')]}, {paths => [file($cwd, 't', 'data', 'file.cat')]}, {paths => [file($cwd, 't', 'data', 'file.txt')]}], 'got correct results for fofn all';

# delimited
ok $ds = VRPipe::DataSource->get(type => 'delimited',
                                 method => 'grouped_single_column',
                                 source => file(qw(t data datasource.fastqs))->absolute->stringify,
                                 options => {delimiter => "\t",
                                             group_by => 1,
                                             column => 2}), 'could create a delimited datasource';

@results = ();
foreach my $element (@{$ds->elements}) {
    push(@results, $element->result);
}
is_deeply \@results, [{paths => [file($cwd, 't/data/2822_6_1.fastq'), file($cwd, 't/data/2822_6_2.fastq')], group => '2822_6'},
                      {paths => [file($cwd, 't/data/2822_7_1.fastq'), file($cwd, 't/data/2822_7_2.fastq')], group => '2822_7'},
                      {paths => [file($cwd, 't/data/2823_4_1.fastq'), file($cwd, 't/data/2823_4_2.fastq')], group => '2823_4'},
                      {paths => [file($cwd, 't/data/8324_8_1.fastq'), file($cwd, 't/data/8324_8_2.fastq')], group => '8324_8'}], 'got correct results for delimited grouped_single_column';

# delimited all columns
ok $ds = VRPipe::DataSource->get(type => 'delimited',
                                 method => 'all_columns',
                                 source => file(qw(t data datasource.2col))->absolute->stringify,
                                 options => {delimiter => "\t"}), 'could create a delimited datasource';

@results = ();
foreach my $element (@{$ds->elements}) {
    push(@results, $element->result);
}
is_deeply \@results, [{paths => [file($cwd, 't/data/file.txt'), file($cwd, 't/data/file2.txt')]},
                      {paths => [file($cwd, 't/data/file3.txt'), file($cwd, 't/data/file.cat')]}], 'got correct results for delimited all_columns';

# test get_methods and method_options
$ds = $ds->_source_instance;
is_deeply [sort $ds->get_methods], [qw(all all_columns grouped_single_column single_column)], 'get_methods returned all the expected methods';
is_deeply [$ds->method_options('single_column')], [['named', 'delimiter', 1, undef, 'Str'],
                                                   ['named', 'column', 1, undef, 'PositiveInt'],
                                                   ['named', 'column_is_path', 0, 1, 'Bool']], 'method_options showed us what options the single_column method takes';

# sequence_index
ok $ds = VRPipe::DataSource->get(type => 'sequence_index',
                                 method => 'lane_fastqs',
                                 source => file(qw(t data datasource.sequence_index))->absolute->stringify,
                                 options => { local_root_dir => $cwd }), 'could create a sequence_index datasource';

@results = ();
foreach my $element (@{$ds->elements}) {
    push(@results, $element->result);
}
is_deeply \@results, [{paths => [file($cwd, 't/data/2822_6.fastq'), file($cwd, 't/data/2822_6_1.fastq'), file($cwd, 't/data/2822_6_2.fastq')], lane => '2822_6'},
                      {paths => [file($cwd, 't/data/2822_7_1.fastq'), file($cwd, 't/data/2822_7_2.fastq')], lane => '2822_7'},
                      {paths => [file($cwd, 't/data/2823_4_1.fastq'), file($cwd, 't/data/2823_4_2.fastq')], lane => '2823_4'},
                      {paths => [file($cwd, 't/data/8324_8_1.fastq'), file($cwd, 't/data/8324_8_2.fastq')], lane => '8324_8'}], 'got correct results for sequence_index lane_fastqs';
my $vrfile = VRPipe::File->get(path => file($cwd, 't/data/2822_6_1.fastq'));
my $meta = $vrfile->metadata;
is_deeply $meta, { expected_md5 => 'f1826489facca0d0bdf02d9586b493f6',
                   lane => '2822_6',
                   study => 'STUDY01',
                   study_name => 'my study name',
                   center_name => 'SC',
                   sample_id => 'SAMPLEID01',
                   sample => 'SAMPLE01',
                   population => 'POP',
                   platform => 'ILLUMINA',
                   library => 'LIB01',
                   insert_size => 200,
                   withdrawn => 0,
                   reads => 200,
                   bases => 12200,
                   analysis_group => 'low coverage',
                   paired => 1,
                   mate => file($cwd, 't/data/2822_6_2.fastq')->stringify }, 'a VRPipe::File created by source has the correct metadata';

# test a special vrtrack test database; these tests are meant for the author
# only, but will also work for anyone with a working VertRes:: and VRTrack::
# setup
SKIP: {
    my $num_tests = 2;
    skip "author-only tests for a VRTrack datasource", $num_tests unless $ENV{VRPIPE_VRTRACK_TESTDB};
    eval "require VertRes::Utils::VRTrackFactory; require VRTrack::VRTrack; require VRTrack::Lane;";
    skip "VertRes::Utils::VRTrackFactory/VRTrack::VRTrack not loading", $num_tests if $@;
    
    # create the vrtrack db
    unless ($ENV{VRPIPE_VRTRACK_TESTDB_READY}) {
        my %cd = VertRes::Utils::VRTrackFactory->connection_details('rw');
        my @sql = VRTrack::VRTrack->schema();
        open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
        print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
        print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
        print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
        foreach my $sql (@sql) {
            print $mysqlfh $sql;
        }
        close($mysqlfh);
        
        # populate it
        system("update_vrmeta.pl --samples t/data/vrtrack.samples --index t/data/vrtrack.sequence.index --database $ENV{VRPIPE_VRTRACK_TESTDB} > /dev/null 2> /dev/null");
        
        # alter processed on the lanes to enable useful tests
        my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
        my $lanes = $vrtrack->processed_lane_hnames();
        my %expectations;
        for my $i (1..60) {
            my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lanes->[$i - 1]);
            my $name = $lane->hierarchy_name;
            
            if ($i <= 50) {
                $lane->is_processed('import' => 1);
                push(@{$expectations{import}}, $name);
                if ($i > 10) {
                    $lane->is_processed('qc' => 1);
                    push(@{$expectations{qc}}, $name);
                }
                if ($i > 20) {
                    $lane->is_processed('mapped' => 1);
                    push(@{$expectations{mapped}}, $name);
                }
                if ($i > 30) {
                    $lane->is_processed('improved' => 1);
                    push(@{$expectations{improved}}, $name);
                }
                if ($i > 40) {
                    $lane->is_processed('stored' => 1);
                    push(@{$expectations{stored}}, $name);
                }
            }
            elsif ($i > 55) {
                $lane->is_withdrawn(1);
                push(@{$expectations{withdrawn}}, $name);
            }
            
            $lane->update;
        }
    }
    
    ok $ds = VRPipe::DataSource->get(type => 'vrtrack',
                                     method => 'lanes',
                                     source => $ENV{VRPIPE_VRTRACK_TESTDB},
                                     options => {import => 1, mapped => 0}), 'could create a vrtrack datasource';
    my $results = 0;
    foreach my $element (@{$ds->elements}) {
        $results++;
    }
    is $results, 20, 'got correct number of results for vrtrack lanes mapped => 0';
   
 
    ### tests for  _has_changed ###
    ok( ! $ds->_source_instance->_has_changed, 'vrtrack datasource _has_changed gives no change' );
    
    # create a new row that has a later time stamp
    my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
    my $name = 'new_lane';
    sleep(2); # wait two seconds to create new lane so we have a later time stamp.
    my $new_lane = VRTrack::Lane->create($vrtrack, $name);
    $new_lane->is_withdrawn(0);
    $new_lane->library_id(16);
    $new_lane->update;
    ok( $ds->_source_instance->_has_changed, 'vrtrack datasource _has_changed gives change after new lane insertion in test vrtrack db');

    # Go back to unchanged state by deleting this lane. Check we don't have any changes
    $new_lane->delete;
    ok( !$ds->_source_instance->_has_changed, 'vrtrack datasource _has_change gives no change after inserted lane deleted in test vrtrack db');

    # add a file to another lane and check for changes
    my $lane_to_add_file_for = VRTrack::Lane->new_by_name($vrtrack, 'ERR003040');
    my $newfile = $lane_to_add_file_for->add_file('new.fastq');
    ok($ds->_source_instance->_has_changed, 'vrtrack datasource _has_changed gives change after adding a file in test vrtrack db');
    $newfile->delete; # return to original state so that _has_changed can be tested again

    # change some md5 sums in the files
    my $file = VRTrack::File->new_by_hierarchy_name( $vrtrack, 'ERR003038.filt.fastq.gz' );
    # check for changes
    $file->md5('34c009157187c5d9a7e976563ec1bad9');
    $file->update;
    ok($ds->_source_instance->_has_changed, 'datasource _has_changed got change after md5 change in file table in test vrtrack db');
    
   ### lane_fastqs tests    
   ok $ds = VRPipe::DataSource->get(type => 'vrtrack',
                                 method => 'lanes_fastqs',
                                 source => $ENV{VRPIPE_VRTRACK_TESTDB},
                                 options => {import => 1, mapped => 0, local_root_dir => file($cwd,'t')->absolute->stringify, library => [  'g1k-sc-NA19190-YRI-1|SC|SRP000542|NA19190'] } ), 'could create a vrtrack datasource';

    @results = ();
    foreach my $element (@{$ds->elements}) {
      push(@results, $element->result);
    }

    is_deeply $results[0], {paths => [file($cwd, 't/data/NA19190/sequence_read/ERR003199.filt.fastq.gz'), file($cwd, 't/data/NA19190/sequence_read/ERR003199_1.filt.fastq.gz'), file($cwd, 't/data/NA19190/sequence_read/ERR003199_2.filt.fastq.gz')], lane => 'ERR003199'}, 'got correct results for vrtrack lane_fastqs'
}
exit;
