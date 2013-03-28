#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 33;
    use VRPipeTest (required_exe => [qw(samtools bcftools fastqcheck)]);
    
    use_ok('VRPipe::FileType');
}

# txt
ok my $ft = VRPipe::FileType->create('txt', { file => file(qw(t data file.txt)) }), 'could create a txt filetype';
is $ft->type, 'txt', 'txt type is correct';
is $ft->check_type(), 1, 'a txt file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [2, 0, 2], 'num_lines, num_header_lines and num_records work for a txt file';

$ft->file(file(qw(t data file.bam)));
is $ft->check_type(), 0, 'a bam file fails the check';

# bam
ok $ft = VRPipe::FileType->create('bam', { file => file(qw(t data file.bam)) }), 'could create a bam filetype';
is $ft->type, 'bam', 'bam type is correct';
is $ft->check_type(), 1, 'a bam file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [92, 89, 3], 'num_lines, num_header_lines and num_records work for a bam file';

throws_ok { $ft = VRPipe::FileType->create('foo', {}); } qr/Invalid implementation class|perhaps you forgot to load/, 'throws when asked to create an invalid filetype';

ok my $parser = $ft->parser, 'could make a parser object';
ok $parser->does('VRPipe::ParserRole'), 'the parser is really a parser';

is $ft->read_backwards, 0, 'default read backwards is off';
# lsf
$ft = VRPipe::FileType->create('lsf', { file => '/foo/bar/baz.lsf' });
is $ft->read_backwards, 1, 'for lsf, read backwards is on';

# vcf
ok $ft = VRPipe::FileType->create('vcf', { file => file(qw(t data file.vcf.gz))->absolute }), 'could create a vcf filetype';
is $ft->type, 'vcf', 'vcf type is correct';
is $ft->check_type(), 1, 'a vcf file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [24, 19, 5], 'num_lines, num_header_lines and num_records work for a vcf, file';

# fq
ok $ft = VRPipe::FileType->create('fq', { file => file(qw(t data 2822_7_2.fastq))->absolute }), 'could create a fq filetype';
is $ft->type, 'fq', 'fq type is correct';
is $ft->check_type(), 1, 'a fq file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [250, 0, 250], 'num_lines, num_header_lines and num_records work for a fq, file';
ok my $fq = VRPipe::File->create(path => file(qw(t data 2822_7_2.fastq))->absolute), 'could create a fq file';
is_deeply [$fq->lines, $fq->lines(raw => 1)], [250, 1000], 'lines and raw lines as expected for fastq file';

# cram
my $ref = file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify;
ok $ft = VRPipe::FileType->create('cram', { file => file(qw(t data hs_chr20.a.cram))->absolute }), 'could create a cram filetype';
is $ft->type, 'cram', 'cram type is correct';
is $ft->check_type(), 1, 'a cram file passes the check';

# bcf
my $bcf = VRPipe::File->create(path => file(qw(t data file.bcf))->absolute);
ok $ft = VRPipe::FileType->create('bcf', { file => $bcf->path }), 'could create a bcf filetype';
is $ft->type, 'bcf', 'bcf type is correct';
is $ft->check_type(),  1, 'a bcf file passes the check';
is $ft->check_magic(), 1, 'bcf magic is correct';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records, $ft->samples], [15671, 35, 15636, [qw(NA12006 NA12144 NA18970 NA18969 NA18995 NA18946)]], 'num_lines, num_header_lines, num_records and samples methods work for a bcf file';

exit;
