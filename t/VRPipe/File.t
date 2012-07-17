#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 36;
    use VRPipeTest;
}

my $vrobj = VRPipe::Manager->get;
my $tmp_dir = $vrobj->tempdir;

my $input_path = file($tmp_dir, 'input.txt');
open(my $fh, '>', $input_path) or die "Could not write to $input_path\n";
print $fh "line1\nline2\n";
close($fh);

ok my $vrfile = VRPipe::File->create(path => $input_path, type => 'txt', metadata => {foo => 'bar'}), 'created a File using create()';
undef($vrfile);
$vrfile = VRPipe::File->get(id => 1);
$vrfile->add_metadata({baz => 'loman'});
is_deeply [$vrfile->path, $vrfile->e, $vrfile->metadata, $vrfile->basename, $vrfile->type, $vrfile->slurp], [$input_path, 1, {foo => 'bar', baz => 'loman'}, 'input.txt', 'txt', "line1\n", "line2\n"], 'file has the expected fields';
cmp_ok $vrfile->s, '>=', 5, 'file has some size';

ok my $orig_mtime = $vrfile->mtime, 'got mtime';
sleep(2);
system("touch ".$vrfile->path);
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
is $vrsource->resolve->id, $real_fileid, 'the source resolves to the final move destination';

$vrdest4->add_metadata({test2 => 'meta2'});
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
is_deeply [$vrdest6->resolve->id, $vrdest4->resolve->id], [$real_fileid, $real_fileid], 'even a move of a symlink still resolves correctly for both source and dest';

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

done_testing;
exit;