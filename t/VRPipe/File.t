#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 17;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPersistentReal;
}

my $vrobj = VRPipe::Manager->get;
my $tmp_dir = $vrobj->tempdir;

my $input_path = file($tmp_dir, 'input.txt');
open(my $fh, '>', $input_path) or die "Could not write to $input_path\n";
print $fh "input1_line1\ninput2_line2\n";
close($fh);

ok my $vrfile = VRPipe::File->get(path => $input_path, type => 'txt', metadata => {foo => 'bar'}), 'created a File using get()';
undef($vrfile);
$vrfile = VRPipe::File->get(id => 1);
$vrfile->add_metadata({baz => 'loman'});
is_deeply [$vrfile->path, $vrfile->e, $vrfile->metadata, $vrfile->basename, $vrfile->type], [$input_path, 1, {foo => 'bar', baz => 'loman'}, 'input.txt', 'txt'], 'file has the expected fields';
cmp_ok $vrfile->s, '>=', 5, 'file has some size';

is $vrfile->raw_lines, 2, 'raw_lines() worked';

my $output_path = file($tmp_dir, 'output.txt');
ok my $vrofile = VRPipe::File->get(path => $output_path, type => 'txt'), 'created another File using get()';
undef($vrofile);
$vrofile = VRPipe::File->get(id => 2);
is_deeply [$vrofile->path, $vrofile->e, $vrofile->s], [$output_path, 0, 0], 'file2 has the expected fields';
throws_ok { VRPipe::File->get(path => 'output.txt', type => 'txt') } qr/must be absolute/, 'using get() with a relative path causes a throw';

ok my $ofh = $vrofile->openw, 'was able to open a file for writing';
print $ofh "foo\n";
$vrofile->close;
ok my $ifh = $vrofile->openr, 'was able to open a file for reading';
my $line = <$ifh>;
$vrofile->close;
is $line, "foo\n", 'the read and write was successfull';

$output_path = file($tmp_dir, 'output.txt.gz');
ok $vrofile = VRPipe::File->get(path => $output_path, type => 'txt'), 'created another File using get()';
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
is $vrofile->raw_lines, 3, 'raw_lines works correctly on a compressed file';

done_testing;
exit;