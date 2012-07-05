#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 6;
    # this test is Sanger-specific, only the author needs to run it
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB)],
                    required_exe => [qw(bamcheck)]);
    use TestPipelines;
    
    use_ok('VRTrack::Factory');
    use_ok('VRPipe::Steps::bamcheck');
    use_ok('VRPipe::Steps::vrtrack_update_mapstats');
}

# setup a little VRTrack db to hold info on our test bams
my %cd = VRTrack::Factory->connection_details('rw');
open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
my @sql = VRTrack::VRTrack->schema();
foreach my $sql (@sql) {
    print $mysqlfh $sql;
}
close($mysqlfh);

my $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
my @lane_names = ('7046_8#8', '7936_7#52');
my %lib_names = ('7046_8#8' => 'KUU25220302_3887679', '7936_7#52' => '5002874');
foreach my $lane_name (@lane_names) {
    my $lane = VRTrack::Lane->create($vrtrack, $lane_name);
    $lane->hierarchy_name($lane_name);
    my $lib = VRTrack::Library->create($vrtrack, $lib_names{$lane_name});
    $lane->library_id($lib->id);
    $lane->update;
}

# vrtrack_auto_qc step requires input files have lane metadata; it also requires
# that it have a populated mapstats table, which in turn means we need the
# bamcheck results stored as metadata on the bams
my $output_dir = get_output_dir('auto_qc');
my $lni = 0;
foreach my $file_prefix (qw(autoqc_normal autoqc_short)) {
    foreach my $suffix ('.bam', '.bam.bamcheck') {
        my $file = VRPipe::File->get(path => file('t', 'data', $file_prefix.$suffix)->absolute);
        $file->add_metadata({lane => $lane_names[$lni]});
    }
    
    # bamcheck metadata
    my $bam = file('t', 'data', $file_prefix.'.bam')->absolute;
    my $bamcheck = file($output_dir, $file_prefix.'.bam.bamcheck');
    VRPipe::Steps::bamcheck->stats_from_bamcheck("bamcheck $bam > $bamcheck");
    
    # mapstats
    VRPipe::Steps::vrtrack_update_mapstats->update_mapstats(db => $ENV{VRPIPE_VRTRACK_TESTDB}, bam => $bam, lane => $lane_names[$lni], plot_dir => $output_dir, plots => []);
    
    $lni++;
}

# autoqc that writes results to vrtrack
my $auto_qc_ps = VRPipe::PipelineSetup->get(name => 'auto_qc',
                                            datasource => VRPipe::DataSource->get(type => 'delimited',
                                                                                  method => 'all_columns',
                                                                                  source => file(qw(t data autoqc.tab))->absolute->stringify,
                                                                                  options => { delimiter => "\t" }),
                                            output_root => $output_dir,
                                            pipeline => VRPipe::Pipeline->get(name => 'vrtrack_auto_qc'),
                                            options => { vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB} });

ok handle_pipeline(), 'vrtrack_auto_qc pipeline ran';

my %actual_auto_qc_files;
$lni = 0;
foreach my $de (@{$auto_qc_ps->datasource->elements}) {
    my @output_subdirs = output_subdirs($de->id, $auto_qc_ps->id);
    my $lane_name = $lane_names[$lni++];
    my $file_name = $lane_name;
    $file_name =~ s/#/_/;
    my $aqcfile = VRPipe::File->get(path => file(@output_subdirs, '1_vrtrack_auto_qc', $file_name.'.auto_qc.txt'));
    if ($aqcfile->s) {
        $actual_auto_qc_files{$lane_name} = [$aqcfile->slurp];
    }
}
my %expected_auto_qc_files;
foreach my $lane_name (@lane_names) {
    my $file_name = $lane_name;
    $file_name =~ s/#/_/;
    my $aqcfile = VRPipe::File->get(path => file('t', 'data', $file_name.'.auto_qc.txt')->absolute);
    $expected_auto_qc_files{$lane_name} = [$aqcfile->slurp];
}
is_deeply \%actual_auto_qc_files, \%expected_auto_qc_files, 'auto qc pipeline generated the expected report txt files showing why the lanes passed';

my $passed_auto_qc_lanes = 0;
my $failed_auto_qc_libs = 0;
$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
foreach my $lane_name (@lane_names) {
    my $lane = VRTrack::Lane->new_by_name($vrtrack, $lane_name);
    $passed_auto_qc_lanes++ if $lane->auto_qc_status eq 'passed';
    my $lib = VRTrack::Library->new($vrtrack, $lane->library_id);
    $failed_auto_qc_libs++ if $lib->auto_qc_status eq 'failed';
}
is_deeply [$passed_auto_qc_lanes, $failed_auto_qc_libs], [0, 2], 'auto qc pipeline set lanes and libraries in VRTrack to passed/failed as appropriate';

finish;