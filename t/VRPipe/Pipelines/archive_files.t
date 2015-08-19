#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Digest::MD5;

BEGIN {
    use Test::Most tests => 10;
    use VRPipeTest (required_env => 'VRPIPE_TEST_PIPELINES');
    use TestPipelines;
}

my $archive_output_dir = get_output_dir('archive_pipeline');
ok my $pipeline = VRPipe::Pipeline->create(name => 'archive_files'), 'able to get the archive_files pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(archive_files);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

# first run a quick pipline to generate some files we can archive
my $test_output_dir = get_output_dir('test_pipeline');
VRPipe::PipelineSetup->create(
    name       => 'my test pipeline setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data datasource.fofn3))->absolute
    ),
    output_root => $test_output_dir,
    pipeline    => VRPipe::Pipeline->create(name => 'test_pipeline'),
    options     => {
        all_option  => 'foo',
        one_option  => 50,
        four_option => 'bar'
    }
);

my @test_files;
my $element_id = 0;
foreach my $in ('file.txt', 'file2.txt', 'file3.txt', 'file4.txt', 'file5.txt', 'file6.txt', 'file7.txt') {
    $element_id++;
    push(@test_files, file(output_subdirs($element_id), '4_test_step_four', "$in.step_one.step_two.step_three.step_four"));
}
ok handle_pipeline(@test_files), 'test pipeline ran and created the output files we will try to archive';

# now run our archive_files pipeline, first creating pool dirs and the pool file
my $dpf = VRPipe::File->create(path => file($archive_output_dir, 'disc_pool_file'), type => 'txt');
my $dpfh = $dpf->openw;
my @pools;
foreach my $pool (qw(pool1 pool2 pool3)) {
    my $dir = dir($archive_output_dir, $pool);
    $pipeline->make_path($dir);
    print $dpfh $dir, "\n";
    push(@pools, $dir);
}
$dpf->close;
my $pool_regex = join('|', @pools);

VRPipe::PipelineSetup->create(
    name       => 'my archive pipeline setup',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'my test pipeline setup[4]',
        options => {}
    ),
    output_root => $archive_output_dir,
    pipeline    => $pipeline,
    options     => { disc_pool_file => $dpf->path->stringify }
);

ok handle_pipeline(), 'archive pipeline ran';
my @archive_files;
foreach my $tfile (@test_files) {
    my $moved_to = VRPipe::File->create(path => $tfile)->resolve->path;
    next unless $moved_to =~ /^$pool_regex/;
    my $expected = file(archive_file_location($tfile))->stringify;
    next unless $moved_to =~ /$expected$/;
    push(@archive_files, $moved_to);
}
is scalar(@archive_files), scalar(@test_files), 'all the files were archived to expected locations';

my $skip_pool;
my @archive_file_paths = map { $_->absolute->stringify } @archive_files;
foreach my $pool (qw(pool1 pool2 pool3)) {
    if (grep /$pool/, @archive_file_paths) {
        $skip_pool = $pool;
        last;
    }
}

# now let's check we can use an altered pool file and it still works
$dpfh  = $dpf->openw;
@pools = ();
foreach my $pool (qw(pool1 pool2 pool3 pool4)) {
    next if $pool eq $skip_pool;
    my $dir = dir($archive_output_dir, $pool);
    $pipeline->make_path($dir);
    print $dpfh $dir, "\n";
    push(@pools, $dir);
}
$dpf->close;
$pool_regex = join('|', @pools);

VRPipe::PipelineSetup->create(
    name       => 'my second archive pipeline setup',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'my test pipeline setup[4]',
        options => {}
    ),
    output_root => $archive_output_dir,
    pipeline    => $pipeline,
    options     => { disc_pool_file => $dpf->path->stringify }
);

ok handle_pipeline(), 'archiving with an altered pool also worked';
my @new_archive_files;
my ($moved_to_path, $moved_from, $md5);
foreach my $tfile (@archive_files) {
    my $moved_to = VRPipe::File->create(path => $tfile)->resolve->path;
    next unless $moved_to =~ /^$pool_regex/;
    
    # only those in $skip_pool should have been archived again
    my $expected = $tfile =~ /$skip_pool/ ? file(archive_file_location($tfile))->stringify : $tfile;
    next unless $moved_to =~ /$expected$/;
    push(@new_archive_files, $moved_to);
    
    if ($tfile =~ /$skip_pool/ && !$moved_from) {
        $moved_to_path = $moved_to;
        $moved_from    = $tfile;
        $md5           = VRPipe::File->get(path => $moved_to)->md5;
    }
}
is scalar(@new_archive_files), scalar(@archive_files), 'all the round-2 files were archived to expected locations';

my $file = VRPipe::File->create(path => $moved_to_path);
$moved_from = VRPipe::File->create(path => $moved_from);
is_deeply [$moved_from->moved_to->id, $file->md5], [$file->id, $md5 ? $md5 : -1], 'moved file has appropriate properties';

VRPipe::PipelineSetup->create(
    name       => 'my third archive pipeline setup',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'my second archive pipeline setup[0]',
        options => {}
    ),
    output_root => $archive_output_dir,
    pipeline    => $pipeline,
    options     => {
        disc_pool_file       => $dpf->path->stringify,
        subdir_from_metadata => 'project_subdir',
    }
);

ok handle_pipeline(), 'archiving to project_subdir worked';
my @subdir_archive_files;
foreach my $tfile (@archive_files) {
    my $moved_to = VRPipe::File->create(path => $tfile)->resolve->path;
    next unless $moved_to =~ /project_subdir/;
    push(@subdir_archive_files, $moved_to);
}
is scalar(@subdir_archive_files), scalar(@archive_files), 'all the round-3 files were archived to expected locations';

finish;
exit;

sub archive_file_location {
    my $file     = VRPipe::File->create(path => shift);
    my $dmd5     = Digest::MD5->new();
    my $ifile_id = $file->id;
    $dmd5->add($ifile_id);
    my $md5 = $dmd5->hexdigest;
    my @chars = split("", $md5);
    return (@chars[0 .. 3], $ifile_id, $file->basename);
}
