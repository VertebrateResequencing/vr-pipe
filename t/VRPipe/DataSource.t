#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class qw(file);
use Cwd;

BEGIN {
    use Test::Most tests => 8;
    
    use_ok('VRPipe::DataSourceFactory');
    
    use TestPersistentReal;
}

ok my $ds = VRPipe::DataSourceFactory->create('list', {method => 'all',
                                                       source => file(qw(t data datasource.list)),
                                                       options => {}}), 'could create a list datasource';

my @results;
while (my $result = $ds->next_result) {
    push(@results, $result);
}
is_deeply \@results, [qw(foo bar henry)], 'got all results correctly';

$ds = VRPipe::DataSourceFactory->create('list', {method => 'all',
                                                 source => file(qw(t data datasource.list)),
                                                 options => {skip_comments => 0}});

@results = ();
while (my $result = $ds->next_result) {
    push(@results, $result);
}
is_deeply \@results, ['foo', 'bar', '# comment', 'henry'], 'got even more results with extra options';

ok $ds = VRPipe::DataSourceFactory->create('fofn', {method => 'all',
                                                    source => file(qw(t data datasource.fofn)),
                                                    options => {}}), 'could create a fofn datasource';

@results = ();
while (my $result = $ds->next_result) {
    push(@results, $result);
}
my $cwd = cwd();
is_deeply \@results, [file($cwd, 't', 'data', 'file.bam'), file($cwd, 't', 'data', 'file.cat'), file($cwd, 't', 'data', 'file.txt')], 'got correct results for fofn all';

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
    
    ok $ds = VRPipe::DataSourceFactory->create('vrtrack', {method => 'lanes',
                                                           source => $ENV{VRPIPE_VRTRACK_TESTDB},
                                                           options => {import => 1, mapped => 0}}), 'could create a vrtrack datasource';
    
    my $results = 0;
    while (my $result = $ds->next_result) {
        $results++;
    }
    is $results, 20, 'got correct number of results for vrtrack lanes mapped => 0';
}

exit;