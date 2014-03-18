#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Spec;
use Parallel::ForkManager;

BEGIN {
    use Test::Most tests => 70;
    use VRPipeTest;
}

my $vrobj      = VRPipe::Manager->get;
my $tmp_dir    = $vrobj->tempdir;
my $vrp_config = VRPipe::Config->new();
my $log_dir    = $vrp_config->testing_logging_directory();
like $tmp_dir, qr/^$log_dir/, 'tempdir defaults to be inside logging dir';
my $fstd = File::Spec->tmpdir;
my $true_tmp_dir = $vrobj->tempdir(DIR => $fstd);
like $true_tmp_dir, qr/^$fstd/, 'tempdir root can be overridden';

my $input_path = file($tmp_dir, 'input.txt');
open(my $fh, '>', $input_path) or die "Could not write to $input_path\n";
print $fh "line1\nline2\n";
close($fh);

ok my $vrfile = VRPipe::File->create(path => $input_path, type => 'txt', metadata => { foo => 'bar' }), 'created a File using create()';
undef($vrfile);
$vrfile = VRPipe::File->get(id => 1);
$vrfile->add_metadata({ baz => 'loman', multi => [qw(val1 val2 val3)] });
is_deeply [$vrfile->path, $vrfile->e, $vrfile->metadata, $vrfile->basename, $vrfile->type, $vrfile->slurp], [$input_path, 1, { foo => 'bar', baz => 'loman', multi => [qw(val1 val2 val3)] }, 'input.txt', 'txt', "line1\n", "line2\n"], 'file has the expected fields';
cmp_ok $vrfile->s, '>=', 5, 'file has some size';
is $vrfile->meta_value('foo'), 'bar', 'meta_value() works';

ok my $orig_mtime = $vrfile->mtime, 'got mtime';
sleep(2);
system("touch " . $vrfile->path);
my $new_mtime = $vrfile->mtime;
is $orig_mtime, $new_mtime, 'mtime unchanged in db';
$vrfile->update_stats_from_disc;
$new_mtime = $vrfile->mtime;
isnt $orig_mtime, $new_mtime, 'mtime updated in db';

is $vrfile->lines, 2, 'lines() worked';

my $output_path = file($tmp_dir, 'output.txt');
ok my $vrofile = VRPipe::File->create(path => $output_path, type => 'txt'), 'created another File using create()';
undef($vrofile);
$vrofile = VRPipe::File->get(id => 2);
is_deeply [$vrofile->path, $vrofile->e, $vrofile->s], [$output_path, 0, 0], 'file2 has the expected fields';
throws_ok { VRPipe::File->create(path => 'output.txt', type => 'txt') } qr/must be absolute/, 'using create() with a relative path causes a throw';

ok my $ofh = $vrofile->openw, 'was able to open a file for writing';
print $ofh "foo\n";
$vrofile->close;
ok my $ifh = $vrofile->openr, 'was able to open a file for reading';
my $line = <$ifh>;
$vrofile->close;
is $line, "foo\n", 'the read and write was successfull';
is $vrofile->lines, 1, 'started with 1 line';
$ofh = $vrofile->open('>>');
print $ofh "2\n3\n4\n";
$vrofile->close;
is_deeply [$vrofile->lines, $vrofile->_lines], [4, 4], 'now has 4 lines after appending';
$ofh = $vrofile->openw;
$vrofile->close;
is_deeply [$vrofile->lines, $vrofile->_lines], [0, undef], 'now has 0 lines after truncate write';

$output_path = file($tmp_dir, 'output.txt.gz');
ok $vrofile = VRPipe::File->create(path => $output_path, type => 'txt'), 'created another File using create()';
ok $ofh = $vrofile->openw, 'was able to open a file for writing compressed';
print $ofh "foo\nbar\nbaz\n";
$vrofile->close;
open(my $gunzip, "gunzip -c $output_path |");
my @lines = <$gunzip>;
close($gunzip);
is_deeply \@lines, ["foo\n", "bar\n", "baz\n"], "we really did write a compressed file";
ok $ifh = $vrofile->openr, 'was able to open a compressed file for reading';
@lines = ();
@lines = <$ifh>;
$vrofile->close;
is_deeply \@lines, ["foo\n", "bar\n", "baz\n"], 'the compressed file could be read with openr';
is $vrofile->lines, 3, 'raw_lines works correctly on a compressed file';

# test copying, moving and symlinking files
my $source_path = file($tmp_dir, 'source.txt');
my $vrsource = VRPipe::File->create(path => $source_path, type => 'txt', metadata => { test => 'meta' });
$ofh = $vrsource->openw;
print $ofh "foo\n";
$vrsource->close;

my $vrcopy = VRPipe::File->create(path => file($tmp_dir, 'copy.txt'));
$vrsource->copy($vrcopy);
is_deeply [$vrsource->e, $vrcopy->e, $vrsource->md5, $vrcopy->md5], [1, 1, 'd3b07384d113edec49eaa6238ad5ff00', 'd3b07384d113edec49eaa6238ad5ff00'], 'copy of a file worked, and it updated the md5s';

my $vrdest1 = VRPipe::File->create(path => file($tmp_dir, 'dest.txt'));
my $vrdest2 = VRPipe::File->create(path => file($tmp_dir, 'dest2.txt'));
my $vrdest3 = VRPipe::File->create(path => file($tmp_dir, 'dest3.txt'));
my $vrdest4 = VRPipe::File->create(path => file($tmp_dir, 'dest4.txt'));
$vrsource->move($vrdest1);
$vrdest1->move($vrdest2);
$vrdest2->symlink($vrdest3);
$vrdest3->symlink($vrdest4);
is_deeply [$vrsource->e, $vrdest1->e, $vrdest2->e, $vrdest3->e, $vrdest4->e, $vrdest4->s], [0, 0, 1, 1, 1, $vrdest2->s], 'file existance and sizes are correct after moves';
is_deeply [$vrdest2->metadata->{test}, $vrdest4->metadata->{test}], ['meta', 'meta'], 'both final moved file and symlink have source metadata';
my $real_fileid = $vrdest2->id;
is $vrdest4->resolve->id, $real_fileid, 'the symlink resolves to the real file';
is $vrdest4->resolve(not_symlinks => 1)->id, $vrdest4->id, 'the symlink resolves to itself in not_symlinks mode';
is $vrsource->resolve->id, $real_fileid, 'the source resolves to the final move destination';

$vrdest4->add_metadata({ test2 => 'meta2' });
$vrdest2->reselect_values_from_db;
is $vrdest2->metadata->{test2}, 'meta2', 'changing metadata on a symlink changes the source file as well';
$ofh = $vrdest4->open('>>');
print $ofh "bar\n";
$vrdest4->close;
$vrdest2->reselect_values_from_db;
is $vrdest2->s, $vrdest4->s, 'writing to a symlink updates the source size';
$vrdest4->update_md5;
$vrdest4->lines;
$vrdest2->reselect_values_from_db;
is_deeply [$vrdest2->md5, $vrdest2->lines], [$vrdest4->md5, $vrdest4->lines], 'updating symlink md5 and lines updates the source md5 and lines as well';

my $vrdest5 = VRPipe::File->create(path => file($tmp_dir, 'dest5.txt'));
$vrdest2->symlink($vrdest5);
$vrdest5->remove;
$vrdest2->reselect_values_from_db;
is_deeply [$vrdest2->md5, $vrdest2->lines], ['f47c75614087a8dd938ba4acff252494', '2'], 'removing a symlink does not undef the md5 or lines of the original file';

my $vrdest6 = VRPipe::File->create(path => file($tmp_dir, 'dest6.txt'));
$vrdest4->move($vrdest6);
ok -l $vrdest6->path, 'moved symlink is a symlink';
is_deeply [$vrdest6->resolve->id, $vrdest4->resolve->id], [$real_fileid, $real_fileid], 'even a move of a symlink still resolves correctly for both source and dest';
my @orig_fids = $vrdest6->original;
is_deeply [\@orig_fids, $vrdest6->original->id], [[$vrdest4->id], $vrdest4->id], 'original() works on moved symlinks';
@orig_fids = $vrdest2->original;
is_deeply [\@orig_fids, $vrdest2->original->id], [[$vrdest1->id, $vrsource->id], $vrsource->id], 'original() works on multiply moved files';

my $vrdest7 = VRPipe::File->create(path => file($tmp_dir, 'dest7.txt'));
my $vrdest8 = VRPipe::File->create(path => file($tmp_dir, 'dest8.txt'));
$vrdest2->symlink($vrdest7);
unlink $vrdest7->path;
$vrdest2->move($vrdest8);
$real_fileid = $vrdest8->id;
$vrdest6->reselect_values_from_db;
my $vrdest6_path = $vrdest6->path;
while (-l $vrdest6_path) { $vrdest6_path = readlink $vrdest6_path; }
is_deeply [$vrdest6_path, $vrdest6->resolve->id], [$vrdest8->path, $real_fileid], 'symlinks were updated to point the the new location of a moved file';
$vrdest7->reselect_values_from_db;
is_deeply [$vrdest7->e, $vrdest7->parent], [0, undef], 'symlink deleted outside of vrpipe was not recreated and existence was updated in db';

my $source_id = $vrsource->id;
undef $vrsource;
$vrfile = VRPipe::File->get(path => $source_path);
is $vrfile->id, $source_id, 'without auto-resolve, getting a VRPipe::File with source path gives you the object for that literal path';
$vrfile = VRPipe::File->get(path => $source_path, auto_resolve => 1);
is $vrfile->id, $real_fileid, 'with auto-resolve, getting a VRPipe::File with source path gives you the object for the moved path';

# test that s() and open() works correctly for relative symlinks
my $symlink = file(qw(t data dirs for symlink higher.link))->absolute;
$vrfile = VRPipe::File->create(path => $symlink);
ok $vrfile->s, 's() worked for a complex relative symlink';
is $vrfile->slurp, "the real file\n", 'slurp also worked on it';

# test that copy/move fail when the destination disk lacks space
$vrsource = VRPipe::File->get(path => $source_path);
my $true_s = $vrsource->s;
$vrsource->s(922337203685477);
$vrsource->update;
$vrcopy->unlink;
throws_ok { $vrsource->copy($vrcopy) } qr/There is not enough disk space available/, 'copy with a file too large for destination throws';
$vrsource->s($true_s);
$vrsource->update;
throws_ok { $vrsource->check_destination_space($vrcopy->path->absolute->dir, 100) } qr/There is not enough disk space remaining/, 'check for there being 100% remaining space at destination throws';

# create_fofn
my $fofn = VRPipe::File->create(path => file($tmp_dir, 'list.fofn'));
$fofn->create_fofn([$vrdest1, $vrdest2, $vrdest3, $vrdest4]);
$fofn->reselect_values_from_db;
is $fofn->lines, 4, 'fofn file created okay';

# test FileMethod's concatenate
my $concat_marker  = "---------------------------------VRPipe--concat---------------------------------\n";
my $skipped_marker = "--[lines skipped during VRPipe-concat]--\n";
$vrsource = VRPipe::File->create(path => file(qw(t data concat.source))->absolute);
my $vrdest = VRPipe::File->create(path => file($tmp_dir, 'concat.dest'));
$vrobj->concatenate($vrsource, $vrdest, unlink_source => 0, add_marker => 1);
is_deeply [$vrdest->slurp], ["line 1\n", "line 2\n", "line 3\n", "line 4\n", "line 5\n", "line 6\n", "line 7\n", "line 8\n", "line 9\n", "line 10\n", $concat_marker], 'concatenate with no limit worked, effectively just copying the file';
$vrdest->unlink;
$vrobj->concatenate($vrsource, $vrdest, unlink_source => 0, add_marker => 1, max_lines => 10);
is_deeply [$vrdest->slurp], ["line 1\n", "line 2\n", "line 3\n", "line 4\n", "line 5\n", "line 6\n", "line 7\n", "line 8\n", "line 9\n", "line 10\n", $concat_marker], 'concatenate with limit of 10';
$vrdest->unlink;
$vrobj->concatenate($vrsource, $vrdest, unlink_source => 0, add_marker => 1, max_lines => 9);
is_deeply [$vrdest->slurp], ["line 1\n", "line 2\n", "line 3\n", "line 4\n", $skipped_marker, "line 6\n", "line 7\n", "line 8\n", "line 9\n", "line 10\n", $concat_marker], 'concatenate with limit of 9';
$vrdest->unlink;
$vrobj->concatenate($vrsource, $vrdest, unlink_source => 0, add_marker => 1, max_lines => 8);
is_deeply [$vrdest->slurp], ["line 1\n", "line 2\n", "line 3\n", "line 4\n", $skipped_marker, "line 7\n", "line 8\n", "line 9\n", "line 10\n", $concat_marker], 'concatenate with limit of 8';

# output_by; we need stepstates and stepoutputfiles to test this, which in turn
# need a bunch of stuff - create all the other object first
my $ds = VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data datasource.fofn3))->absolute);
$ds->elements;
my $pipeline = VRPipe::Pipeline->create(name => 'archive_files');
my $setup = VRPipe::PipelineSetup->create(
    name        => 'my archive pipeline setup',
    datasource  => $ds,
    output_root => '/tmp',
    pipeline    => $pipeline,
    options     => { disc_pool_file => '/tmp' },
    active      => 0
);
my $stepstate  = VRPipe::StepState->create(stepmember => 1, dataelement => 1, pipelinesetup => 1, complete => 1);
my $stepstate2 = VRPipe::StepState->create(stepmember => 1, dataelement => 2, pipelinesetup => 1, complete => 1);
my $step_out_file  = VRPipe::File->create(path => "/step/o/file");
my $step_out_file2 = VRPipe::File->create(path => "/step/o/file2");
my $sof  = VRPipe::StepOutputFile->create(stepstate => $stepstate->id,  file => $step_out_file->id,  output_key => "foo");
my $sof2 = VRPipe::StepOutputFile->create(stepstate => $stepstate2->id, file => $step_out_file2->id, output_key => "foo");
# archive_files pipeline moves files without specifying the result was an
# output file; test that we at least know the result was an output of the
# original setup that archive_files was run on
my $moved_out_file = VRPipe::File->create(path => "/moved/o/file");
$step_out_file->moved_to($moved_out_file->id);
$step_out_file->update;
my $scalar_context_result = $moved_out_file->output_by();
is $scalar_context_result, 1, 'output_by in scalar context returns true for a moved output file';
my @list_context_result = $moved_out_file->output_by();
is_deeply [map { $_->id } @list_context_result], [$stepstate->id], 'output_by in list context returns the correct stepstate';
# test that we can see when multiple stepstates created the same file
my $stepstate3 = VRPipe::StepState->create(stepmember => 1, dataelement => 3, pipelinesetup => 1, complete => 1, same_submissions_as => $stepstate->id);
my $sof3 = VRPipe::StepOutputFile->create(stepstate => $stepstate3->id, file => $step_out_file->id, output_key => "foo");
@list_context_result = $step_out_file->output_by();
is_deeply [sort map { $_->id } @list_context_result], [$stepstate->id, $stepstate3->id], 'we see both stepstates when both created the same file';
my $special = $step_out_file->output_by(1);
is $special->id, $stepstate->id, 'and in the special mode we see the first, since the other had the same subs';
# they might not always have same_submissions_as, so test when they do and do
# not have the same job
my $stepstate4 = VRPipe::StepState->create(stepmember => 1, dataelement => 4, pipelinesetup => 1, complete => 1);
my $sof4 = VRPipe::StepOutputFile->create(stepstate => $stepstate4->id, file => $step_out_file->id, output_key => "foo");
my $job1 = VRPipe::Job->create(cmd => 'foo');
my $job2 = VRPipe::Job->create(cmd => 'bar');
my $req = VRPipe::Requirements->create(memory => 1, time => 1);
my $sub1 = VRPipe::Submission->create(job => $job1->id, stepstate => $stepstate->id,  requirements => $req->id);
my $sub2 = VRPipe::Submission->create(job => $job2->id, stepstate => $stepstate4->id, requirements => $req->id);
@list_context_result = $step_out_file->output_by();
is_deeply [sort map { $_->id } @list_context_result], [$stepstate->id, $stepstate3->id, $stepstate4->id], 'we see all 3 stepstates when we added another';
$special = $step_out_file->output_by(1);
is $special, undef, 'but in the special mode we see none of them, since one of them had a different job';
$sub2->job($job1->id);
$sub2->update;
$special = $step_out_file->output_by(1);
is $special->id, $stepstate->id, 'and in the special mode, when it has no same_submissions_as but does have the same job, we see the first stepstate again';
# test it works when the output file was created as a symlink and then moved
my $ifile = VRPipe::File->get(path => file(qw(t data file.txt))->absolute);
my $step_out_symlink = VRPipe::File->create(path => "/step/o/symlink", parent => $ifile->id);
my $stepstate5 = VRPipe::StepState->create(stepmember => 1, dataelement => 5, pipelinesetup => 1, complete => 1);
my $sof5 = VRPipe::StepOutputFile->create(stepstate => $stepstate5->id, file => $step_out_symlink->id, output_key => "foo");
is_deeply [map { $_->id } $step_out_symlink->output_by], [$stepstate5->id], 'output_by works on a symlink output';
my $moved_out_symlink = VRPipe::File->create(path => '/moved/o/symlink', parent => $ifile->id);
$step_out_symlink->moved_to($moved_out_symlink->id);
$step_out_symlink->update;
is_deeply [map { $_->id } $moved_out_symlink->output_by], [$stepstate5->id], 'output_by works on a symlink output that was moved';

# input_to
my $ifile7 = VRPipe::File->get(path => file(qw(t data file7.txt))->absolute);
is $step_out_file->input_to, 0, 'input_to() on a non dataelement file returns 0';
is $ifile->input_to,         1, 'input_to() on a dataelement file returns 1';
is_deeply [[map { $_->id } $ifile->input_to], [map { $_->id } $ifile7->input_to]], [[1], [7]], 'input_to works in list context';
my $moved_ifile = VRPipe::File->create(path => '/tmp/moved_ifile.txt');
$ifile->moved_to($moved_ifile->id);
$ifile->update;
is_deeply [map { $_->id } $moved_ifile->input_to], [1], 'input_to works on a moved input file';

# test that we can add different metadata in 10 roughly overlapping processes (2
# different updates to 5 different files) and all metadata is kept
my $fm = Parallel::ForkManager->new(10);
my $common_meta = { original_pg_chain => 'asdlfk asdlfkj a;lskfj ;asdkfj;laksd fj;aksj dfl;kjasdf ;laskdjf;lkaj sdf;' };
my %unique_meta;
my @meta_files;
for my $i (1 .. 5) {
    for my $l (qw(a b c d e f g h i j k)) {
        $unique_meta{$i}->{$l} = $l . $i;
    }
    push(@meta_files, VRPipe::File->create(path => file($tmp_dir, 'metafile.' . $i)));
}
for my $i (1 .. 5) {
    my $meta_file = $meta_files[$i - 1];
    
    for my $j (1 .. 2) {
        my $meta = $j == 1 ? $common_meta : $unique_meta{$i};
        
        $fm->start and next;
        
        sleep($j); # this is critical for revealing the bug
        $meta_file->add_metadata($meta);
        
        $fm->finish(0);
    }
}
$fm->wait_all_children;

my %actual_meta;
my %expected_meta;
for my $i (1 .. 5) {
    my $meta_file = $meta_files[$i - 1];
    $meta_file->reselect_values_from_db;
    $actual_meta{$i} = $meta_file->metadata;
    $expected_meta{$i} = { %{ $unique_meta{$i} }, %{$common_meta} };
}

is_deeply \%actual_meta, \%expected_meta, 'add_metadata() worked for 10 processes roughly overlapping';

my @kvlms = VRPipe::KeyValListMember->search({ keyval_key => 'original_pg_chain' });
is scalar(@kvlms), 6, 'No KeyValListMembers were duplicated incorrectly';

# also test that we don't duplicate identical filelists when they're created
# in parallel
my @fillist_files;
for my $i (1 .. 10) {
    push(@fillist_files, VRPipe::File->create(path => file($tmp_dir, 'filelist.' . $i)));
}
for my $i (1 .. 10) {
    $fm->start and next;
    
    VRPipe::FileList->create(files => \@fillist_files);
    
    $fm->finish(0);
}
$fm->wait_all_children;

my @fls = VRPipe::FileList->search({});
is scalar(@fls), 8, 'only 1 FileList was created when creating the same list 10 times in parallel'; # the other 7 were made by previous tests

# test the common_metadata method in StepRole; create a Step to get access to
# the method
my $step = VRPipe::Step->create(
    name               => 'step',
    inputs_definition  => {},
    body_sub           => sub { return 1; },
    outputs_definition => {},
    post_process_sub   => sub { return 1 },
    description        => 'step'
);
my %common = (foo => 'bar', cat => 'dog', multi => [qw(ay be ce)]);
@meta_files = ();
for my $i (1 .. 5) {
    push(@meta_files, VRPipe::File->create(path => file($tmp_dir, 'metafileb.' . $i), metadata => { %common, unique => $i }));
}
is_deeply $step->common_metadata(\@meta_files), \%common, 'common_metadata() from StepRole works, excluding a unique scalar';
$meta_files[0]->add_metadata({ multi => [qw(ay be ce de)] });
is_deeply $step->common_metadata(\@meta_files), { foo => 'bar', cat => 'dog' }, 'common_metadata() from StepRole works, excluding a multi value key that differs';

done_testing;
exit;
