#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class qw(file);

BEGIN {
    use Test::Most tests => 46;
    use VRPipeTest (required_env => 'SAMTOOLS');
    
    use_ok('VRPipe::Parser');
}

my $pb = VRPipe::Parser->create('bam', { file => file(qw(t data parser.bam)) });

# first test all the flag-related methods, which are independent of the
# particular bam
my $flag = 163;
is $pb->is_sequencing_paired($flag),   1, 'is_sequencing_paired test';
is $pb->is_mapped_paired($flag),       1, 'is_mapped_paired test';
is $pb->is_mapped($flag),              1, 'is_mapped test';
is $pb->is_mate_mapped($flag),         1, 'is_mate_mapped test';
is $pb->is_reverse_strand($flag),      0, 'is_reverse_strand test';
is $pb->is_mate_reverse_strand($flag), 1, 'is_mate_reverse_strand test';
is $pb->is_first($flag),               0, 'is_first test';
is $pb->is_second($flag),              1, 'is_second test';
is $pb->is_second(0), 0, 'is_second test on flag 0';
is $pb->is_primary($flag),   1, 'is_primary test';
is $pb->passes_qc($flag),    1, 'passes_qc test';
is $pb->is_duplicate($flag), 0, 'is_duplicate test';

# now test parsing actual bam records; header methods are tested in Parser.t
ok my $pr = $pb->parsed_record(), 'parsed_record returned something';
is ref($pr), 'HASH', 'parsed_record returns a hash ref';
is keys %{$pr}, 0, 'the parsed_record starts off empty';

$pb->get_fields('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM', 'MPOS', 'ISIZE', 'SEQ', 'QUAL');
ok $pb->next_record, 'next_record now works';
is_deeply $pr, { QNAME => 'IL3_2822:6:1:3:1989', FLAG => '69', RNAME => '*', POS => '0', MAPQ => '0', CIGAR => '*', MRNM => '*', MPOS => '0', ISIZE => '0', SEQ => 'AAAAAAAAAAAAAAAAAAAAAAAAAAANNAAAAAAAANAANAAAAAAAAAAAAAAAAAAAA', QUAL => '44444477774774774447477747($$(447474\'$(($(7777477744477774477' }, 'parsed_record contains correct info for first line';
$pb->get_fields('SEQ');
ok $pb->next_record, 'next_record worked again';
is $pr->{SEQ}, 'NANANNAAANNANNANAAAAAAAAANAAANNNNNNNNNNNNNNNNNNANNNNNN', 'parsed_record contains correct info for second line';
is $pr->{QNAME}, undef, 'parsed_record contains undef for an unrequested field';

# check the last line as well
$pb->get_fields('QNAME', 'NM', 'MD');
while ($pb->next_record) {
    next;
}
is $pr->{QNAME}, 'IL3_2822:6:1:46:1716', 'parsed_record contains correct qname for last line';
is $pr->{NM},    1,                      'parsed_record contains correct NM for last line';
is $pr->{MD},    '47A6',                 'parsed_record contains correct MD for last line';

# more tests in a different file
my $b_file = Path::Class::File->new(qw(t data parser.2.bam));
$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->get_fields('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM', 'MPOS', 'ISIZE', 'SEQ', 'QUAL', 'XT', 'NM', 'SM', 'AM', 'X0', 'X1', 'XM', 'XO', 'XG', 'MD', 'SEQ_LENGTH', 'MAPPED_SEQ_LENGTH');
$pb->next_record;
$pr = $pb->parsed_record();
is_deeply $pr, { QNAME => 'IL3_2822:6:1:40:1780', FLAG => 163, RNAME => 'Streptococcus_suis', POS => 4927, MAPQ => 37, CIGAR => '54M', MRNM => 'Streptococcus_suis', MPOS => 5086, ISIZE => 213, SEQ => 'CTAGAAGACGGAAAATCTGCCCGCACAGTTGAGTTCACAGATGAAGAACAAAAA', QUAL => '?>??;:@<?????6)5><?8=;??=;;?2:?>>2<>:>:?:3@97?1764091=', XT => 'U', NM => 0, SM => 37, AM => 37, X0 => 1, X1 => 0, XM => 0, XO => 0, XG => 0, MD => 54, SEQ_LENGTH => 54, MAPPED_SEQ_LENGTH => 54 }, 'next_record test on first line';
$pb->next_record;
is $pr->{ISIZE}, -213, 'get_fields second line mate gets a negative isize';
$pb->get_fields('QNAME', 'FLAG', 'RG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM', 'MPOS', 'ISIZE', 'SEQ', 'SEQ_LENGTH', 'QUAL', 'XT');

while ($pb->next_record) {
    next;
}
is_deeply $pr, { QNAME => 'IL3_2822:6:1:46:1362', FLAG => 133, RG => '*', RNAME => '*', POS => 0, MAPQ => 0, CIGAR => '*', MRNM => '*', MPOS => 0, ISIZE => 0, SEQ => 'CGTTGTAGCTGAGGCTGACATTGAATTTTCACATCTCGAAAGCTTCATCAGTCT', SEQ_LENGTH => 54, QUAL => '??,<-(\'4+9/&<A>8((2)4<?835?>.9+(.\'%39@<8<3\'6,.)-==./.(', XT => '*' }, 'next_record test on last line';

# test region();
$pb = VRPipe::Parser->create('bam', { file => $b_file });
my $c = 0;
foreach my $region ('Streptococcus_suis:29888-50000', 'Streptococcus_suis:80000-90000', 'Streptococcus_suis:60000-70000') {
    $pb->region($region);
    while ($pb->next_record) {
        $c++;
    }
}
is $c, 18, 'getting multiple regions test';

# test writing
my $for_tmp_dir = VRPipe::Manager->get;
my $temp_dir    = $for_tmp_dir->tempdir();
my $out_bam1    = Path::Class::File->new($temp_dir, 'out1.bam');
my $out_bam2    = Path::Class::File->new($temp_dir, 'out2.bam');
my $out_bam3    = Path::Class::File->new($temp_dir, 'out3.bam');
$pb->region('');
while ($pb->next_record) {
    $pb->write_result($out_bam1);
}
$pb->close;
is_deeply [get_bam_header($out_bam1)], [get_bam_header($b_file)], 'having set input region to empty string, header of an output file matches the input file';
is scalar(get_bam_body($out_bam1)), 2000, 'output bam has correct number of lines';

$b_file = Path::Class::File->new(qw(t data parser.1kg_lane.bam));
$pb     = VRPipe::Parser->create('bam', { file => $b_file });
$c      = 0;
$pb->ignore_tags_on_write(qw(OQ XM XG XO));
while ($pb->next_record) {
    $c++;
    if ($c == 1) {
        $pb->write_result($out_bam2);
    }
    elsif ($c == 2) {
        $pb->write_result($out_bam3);
    }
    else {
        last;
    }
}
undef $pb;
my @g1k_lines = grep { s/\tOQ:\S+|\tXM:\S+|\tXG:\S+|\tXO:\S+//g } get_bam_body($b_file);
is_deeply [get_bam_body($out_bam2)], [$g1k_lines[0]], 'ignore_tags_on_write has its effect';
is_deeply [get_bam_body($out_bam3)], [$g1k_lines[1]], 'we can output multiple files in a single next_record loop';

# filtering methods
$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->required_readgroup('SRR035022');
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 3, 'required_readgroup with the correct RG gives all results';
$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->required_readgroup('foo');
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 0, 'required_readgroup with an unused RG gives no results';
$pb = VRPipe::Parser->create('bam', { file => $b_file });
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 3, 'no required_readgroup in a new instance gives us all results';

$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->required_library('foo');
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 0, 'required_library with an unused lib gives no results';
$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->required_library('Solexa-16652');
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 3, 'required_library with the correct lib gives all results';

$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->required_flag(32);
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 2, 'required_flag gives a subset of results';
$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->filtering_flag(32);
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 1, 'filtering_flag gives the opposite subset of results';

$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->flag_selector(mate_reverse => 1, '1st_in_pair' => 1);
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 2, 'flag_selector worked with two true flag names';
$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->flag_selector(mate_reverse => 0, '1st_in_pair' => 0);
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 0, 'flag_selector worked with two false flag names';
$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->flag_selector(mate_reverse => 0, '1st_in_pair' => 1);
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 1, 'flag_selector worked with a mix of true and false flag names';

$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->minimum_quality(30);
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 0, 'minimum_quality gives a subset of results';
$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->minimum_quality(0);
$c = 0;
while ($pb->next_record) {
    $c++;
}
is $c, 3, 'minimum_quality(0) gives all results';

# hard/soft clipping seq length tests
$b_file = Path::Class::File->new(qw(t data parser.hard_soft.bam));
$pb = VRPipe::Parser->create('bam', { file => $b_file });
$pb->get_fields('SEQ', 'SEQ_LENGTH', 'MAPPED_SEQ_LENGTH');
$pr = $pb->parsed_record;
$pb->next_record;
is_deeply $pr, { SEQ => 'CCCATAGCCCTATCCCTAACCCTAACCCGAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC', SEQ_LENGTH => 76, MAPPED_SEQ_LENGTH => 72 }, 'seq length correct both raw and clipped';

# double-check we can parse OQ tags in improved bams
$pb = VRPipe::Parser->create('bam', { file => file(qw(t data 2822_6.improved.pe.bam)) });
$pb->get_fields('QNAME', 'OQ');
$pb->next_record;
$pr = $pb->parsed_record();
is_deeply $pr, { QNAME => 'IL3_2822:6:1:20:1902', OQ => '?>?<?>/445224;97556294/359543443595:79998649;7<=999?>?=B@==>>' }, 'able to get OQ out of an improved bam';

exit;

sub get_bam_header {
    my $bam_file = shift;
    open(my $bamfh, "samtools view -H $bam_file |") || die "Could not open samtools view -H $bam_file\n";
    my @header_lines;
    while (<$bamfh>) {
        chomp;
        push(@header_lines, $_);
    }
    close($bamfh);
    return @header_lines;
}

sub get_bam_body {
    my $bam_file = shift;
    open(my $bamfh, "samtools view $bam_file |") || die "Could not open samtools view $bam_file\n";
    my @records;
    while (<$bamfh>) {
        chomp;
        next if /^@/;
        push(@records, $_);
    }
    close($bamfh);
    return @records;
}
