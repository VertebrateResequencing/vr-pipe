#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 63;
    use VRPipeTest (required_exe => [qw(samtools fastqcheck htsfile)]);
    
    use_ok('VRPipe::FileType');
}

# txt
ok my $ft = VRPipe::FileType->create('txt', { file => file(qw(t data file.txt)) }), 'could create a txt filetype';
isa_ok($ft, 'VRPipe::FileType::txt');
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

ok my $parser = $ft->parser, 'could make a bam parser object';
ok $parser->does('VRPipe::ParserRole'), 'the bam parser is really a parser';

# cram
ok $ft = VRPipe::FileType->create('cram', { file => file(qw(t data hs_chr20.a.cram))->absolute }), 'could create a cram filetype';
is $ft->type, 'cram', 'cram type is correct';
is $ft->check_type(), 1, 'a cram file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [917, 5, 912], 'num_lines, num_header_lines and num_records work for a bam file of type aln';

ok $parser = $ft->parser, 'could make a cram parser object';
ok $parser->does('VRPipe::ParserRole'), 'the cram parser is really a parser';

# aln -- bam or cram
ok $ft = VRPipe::FileType->create('aln', { file => file(qw(t data file.bam)) }), 'could create a aln filetype from bam file';
is $ft->type, 'aln', 'aln type is correct for bam file';
is $ft->check_type(), 1, 'a bam file passes the aln file type check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [92, 89, 3], 'num_lines, num_header_lines and num_records work for a bam file of type aln';

ok $ft = VRPipe::FileType->create('aln', { file => file(qw(t data hs_chr20.a.cram))->absolute }), 'could create a aln filetype from cram file';
is $ft->type, 'aln', 'aln type is correct for cram file';
is $ft->check_type(), 1, 'a cram file passes the aln file type check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [917, 5, 912], 'num_lines, num_header_lines and num_records work for a cram file of type aln';

# check we don't validate just any old thing compressed with bgzip
ok $ft = VRPipe::FileType->create('bam', { file => file(qw(t data bgzip_compressed_text_file.gz))->absolute }), 'could create a bam filetype with a bgzipped file';
is $ft->type, 'bam', 'its type is bam';
is $ft->check_type(), 0, 'the bgzipped file fails to validate as a bam';

is $ft->read_backwards, 0, 'default read backwards is off';
# lsf
$ft = VRPipe::FileType->create('lsf', { file => '/foo/bar/baz.lsf' });
is $ft->read_backwards, 1, 'for lsf, read backwards is on';

# vcf
ok $ft = VRPipe::FileType->create('vcf', { file => file(qw(t data file.vcf.gz))->absolute }), 'could create a vcf filetype';
is $ft->type, 'vcf', 'vcf type is correct';
is $ft->check_type(), 1, 'a vcf file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records, $ft->samples, $ft->num_samples], [25, 20, 5, [qw(NA00001 NA00002 NA00003)], 3], 'num_lines, num_header_lines, num_records, samples and num_samples work for a vcf file';

# bcf
my $bcf = VRPipe::File->create(path => file(qw(t data file.bcf))->absolute);
ok $ft = VRPipe::FileType->create('bcf', { file => $bcf->path }), 'could create a bcf filetype';
is $ft->type, 'bcf', 'bcf type is correct';
is $ft->check_type(), 1, 'a bcf file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records, $ft->samples, $ft->num_samples], [25, 20, 5, [qw(NA00001 NA00002 NA00003)], 3], 'num_lines, num_header_lines, num_records, samples and num_samples work for a bcf file';

# var -- vcf or bcf
ok $ft = VRPipe::FileType->create('var', { file => file(qw(t data file.vcf.gz))->absolute }), 'could create a var filetype from vcf';
is $ft->type, 'var', 'var type is correct for vcf';
is $ft->check_type(), 1, 'a vcf file passes the var filetype check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records, $ft->samples, $ft->num_samples], [25, 20, 5, [qw(NA00001 NA00002 NA00003)], 3], 'num_lines, num_header_lines, num_records, samples and num_samples work for a vcf file of type var';

ok $ft = VRPipe::FileType->create('var', { file => $bcf->path }), 'could create a var filetype from bcf';
is $ft->type, 'var', 'var type is correct for bcf';
is $ft->check_type(), 1, 'a bcf file passes the var filetype check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records, $ft->samples, $ft->num_samples], [25, 20, 5, [qw(NA00001 NA00002 NA00003)], 3], 'num_lines, num_header_lines, num_records, samples and num_samples work for a bcf file of type var';

# fq
ok $ft = VRPipe::FileType->create('fq', { file => file(qw(t data 2822_7_2.fastq))->absolute }), 'could create a fq filetype';
is $ft->type, 'fq', 'fq type is correct';
is $ft->check_type(), 1, 'a fq file passes the check';
is_deeply [$ft->num_lines, $ft->num_header_lines, $ft->num_records], [250, 0, 250], 'num_lines, num_header_lines and num_records work for a fq, file';
ok my $fq = VRPipe::File->create(path => file(qw(t data 2822_7_2.fastq))->absolute), 'could create a fq file';
is_deeply [$fq->lines, $fq->lines(raw => 1)], [250, 1000], 'lines and raw lines as expected for fastq file';

# check we don't validate just any old thing compressed with bgzip
ok $ft = VRPipe::FileType->create('bcf', { file => file(qw(t data bgzip_compressed_text_file.gz))->absolute }), 'could create a bcf filetype with a bgzipped file';
is $ft->type, 'bcf', 'its type is bcf';
is $ft->check_type(), 0, 'the bgzipped file fails to validate as a bcf';

# arbitrary filetypes should be supported and treated like 'any', just checking
# the file extension
my $type1 = VRPipe::File->create(path => file(qw(t data 2822_6_1.typ1))->absolute);
ok $ft = VRPipe::FileType->create('typ1', { file => $type1->path }), 'could create a typ1 filetype';
isa_ok($ft, 'VRPipe::FileType::any');
is $ft->type, 'typ1', 'typ1 type is correct';
is $ft->check_type(), 1, 'a typ1 file passes the check';
$ft->file(file(qw(t data 2822_6.pe.typ2))->absolute);
is $ft->check_type(), 0, 'a typ2 file failes a typ1 check';

# because the database has a length limit of 8, however, we don't allow long
# filetypes
throws_ok { $ft = VRPipe::FileType->create('foobarbaz', {}); } qr/Invalid implementation class|perhaps you forgot to load/, 'throws when asked to create an invalid filetype of greater than 8 characters';

exit;
