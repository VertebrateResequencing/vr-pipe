#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 8;
    
    use_ok('VRPipe::FileType');
}

ok my $ft = VRPipe::FileType->create('txt', {file => '/foo/bar/baz.txt'}), 'could create a txt filetype';
is $ft->type, 'txt', 'type is correct';
is $ft->check_type(), 1, 'a txt file passes the check';
$ft->file('/foo/bar/baz.bam');
is $ft->check_type(), 0, 'a bam file fails the check';

ok $ft = VRPipe::FileType->create('bam', {file => '/foo/bar/baz.bam'}), 'could create a bam filetype';
is $ft->type, 'bam', 'type is correct';

throws_ok {$ft = VRPipe::FileType->create('foo', {});} qr/Invalid implementation class/, 'throws when asked to create an invalid filetype';

exit;