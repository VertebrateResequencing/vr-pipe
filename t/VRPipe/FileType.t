#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 13;
    
    use_ok('VRPipe::FileType');
    
    use TestPersistentReal;
}

ok my $ft = VRPipe::FileType->create('txt', {file => 't/data/file.txt'}), 'could create a txt filetype';
is $ft->type, 'txt', 'type is correct';
is $ft->check_type(), 1, 'a txt file passes the check';
$ft->file('t/data/file.bam');
is $ft->check_type(), 0, 'a bam file fails the check';

ok $ft = VRPipe::FileType->create('bam', {file => 't/data/file.bam'}), 'could create a bam filetype';
is $ft->type, 'bam', 'type is correct';
is $ft->check_type(), 1, 'a bam file passes the check';

throws_ok {$ft = VRPipe::FileType->create('foo', {});} qr/Invalid implementation class/, 'throws when asked to create an invalid filetype';

ok my $parser = $ft->parser, 'could make a parser object';
ok $parser->does('VRPipe::ParserRole'), 'the parser is really a parser';

is $ft->read_backwards, 0, 'default read backwards is off';
$ft = VRPipe::FileType->create('lsf', {file => '/foo/bar/baz.lsf'});
is $ft->read_backwards, 1, 'for lsf, read backwards is on';


exit;