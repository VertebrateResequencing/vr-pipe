#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 22;
    use VRPipeTest (required_exe => [qw(samtools)]);
    
    use_ok('VRPipe::FileType');
}

# txt
ok my $ft = VRPipe::FileType->create('txt', { file => file(qw(t data file.txt)) }), 'could create a txt filetype';
is $ft->type, 'txt', 'type is correct';
is $ft->check_type(), 1, 'a txt file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [2, 0, 2], 'num_lines, num_header_lines and num_records work for a txt file';

$ft->file(file(qw(t data file.bam)));
is $ft->check_type(), 0, 'a bam file fails the check';

# bam
ok $ft = VRPipe::FileType->create('bam', { file => file(qw(t data file.bam)) }), 'could create a bam filetype';
is $ft->type, 'bam', 'type is correct';
is $ft->check_type(), 1, 'a bam file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [92, 89, 3], 'num_lines, num_header_lines and num_records work for a bam, file';

throws_ok { $ft = VRPipe::FileType->create('foo', {}); } qr/Invalid implementation class|perhaps you forgot to load/, 'throws when asked to create an invalid filetype';

ok my $parser = $ft->parser, 'could make a parser object';
ok $parser->does('VRPipe::ParserRole'), 'the parser is really a parser';

is $ft->read_backwards, 0, 'default read backwards is off';
# lsf
$ft = VRPipe::FileType->create('lsf', { file => '/foo/bar/baz.lsf' });
is $ft->read_backwards, 1, 'for lsf, read backwards is on';

# vcf
ok $ft = VRPipe::FileType->create('vcf', { file => file(qw(t data file.vcf.gz))->absolute }), 'could create a vcf filetype';
is $ft->type, 'vcf', 'type is correct';
is $ft->check_type(), 1, 'a vcf file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [24, 19, 5], 'num_lines, num_header_lines and num_records work for a vcf, file';

# cram
my $ref = file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify;
ok $ft = VRPipe::FileType->create('cram', { file => file(qw(t data hs_chr20.a.cram))->absolute }), 'could create a cram filetype';
is $ft->type, 'cram', 'type is correct';
is $ft->check_type(), 1, 'a cram file passes the check';

exit;
