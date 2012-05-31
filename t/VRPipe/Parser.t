#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class qw(file);

BEGIN {
    use Test::Most tests => 74;
    use VRPipeTest;
    
    use_ok('VRPipe::Parser');
}

# lsf
ok my $p = VRPipe::Parser->create('lsf', {file => file(qw(t data lsf.stdout))}), 'could create an lsf parser';
$p->fh;
is_deeply [$p->file, $p->_vrpipe_file->path], [file(qw(t data lsf.stdout)), file(qw(t data lsf.stdout))->absolute], 'file details are correct';

throws_ok {$p = VRPipe::Parser->create('foo', {});} qr/Invalid implementation class|perhaps you forgot to load/, 'throws when asked to create an invalid parser';

is $p->memory, 1, 'memory is correct without a next_record() call';
is $p->cpu_time, 4.75, 'so is cpu_time';

ok $p->next_record, 'next_record worked';
is_deeply $p->parsed_record, [q[perl -MVRPipe::Persistent::SchemaBase -MVRPipe::Scheduler -e "VRPipe::Persistent::SchemaBase->database_deployment(q[testing]); VRPipe::Scheduler->get(id => 3)->run_on_node(index => shift, array => 4);"],
                              'OK',
                              1,
                              13,
                              4.98,
                              0.38,
                              'normal',
                              892564,
                              8], 'parsed_record contains all the correct details';

ok $p->next_record, 'next_record worked again';
is $p->cpu_time, 5.05, 'cpu_time worked after next_record calls';

ok ! $p->next_record, 'next_record returns false when no more records';

undef $p;
$p = VRPipe::Parser->create('lsf', {file => file(qw(t data lsf.stdout))});
$p->next_record;
is $p->parsed_record->[4], 4.75, 'next_record and parsed_record work on the first (last) record';

# loc
undef $p;
$p = VRPipe::Parser->create('loc', {file => file(qw(t data parser.loc))});
is $p->time, 20, 'local parser time is correct';
is $p->cmd, 'perl -MVRPipe::Persistent::Schema -e "VRPipe::Persistent::SchemaBase->database_deployment(q[testing]); VRPipe::Scheduler->get(id => 1)->run_on_node(index => shift, array => 2);"', 'local parser cmd is correct';
is $p->status, 'OK', 'local parser status is correct';
is $p->sid, 2, 'local parser sid is correct';
is $p->aid, 1, 'local parser aid is correct';

throws_ok {$p = VRPipe::Parser->create('cat', {file => 't/data/not.extant'}); $p->next_record; } qr/does not exist, so cannot be opened for reading/, 'throws when asked to parse a non-existant file';

# cat
undef $p;
$p = VRPipe::Parser->create('cat', {file => file(qw(t data file.cat))});
my @records;
while ($p->next_record) {
    push(@records, [@{$p->parsed_record}]);
}
is_deeply [@records], [["first line of 4th record", "second line of 4th record"],
                       [],
                       ["first line of 2nd record", "second line of 2nd record"],
                       ["first line of 1st record", "second line of 1st record"]], 'cat file was parsed correctly';
$p = VRPipe::Parser->create('cat', {file => file(qw(t data parser.cat))});
@records = ();
while ($p->next_record) {
    push(@records, [@{$p->parsed_record}]);
}
is_deeply [@records], [[],
                       ["stdout message: failing on purpose since this is try 2"],
                       ["stdout message: failing on purpose since this is try 1"]], 'another cat file was parsed correctly';

# fqc
$p = VRPipe::Parser->create('fqc', {file => file(qw(t data parser.fastqcheck ))});
is_deeply [$p->num_sequences, $p->total_length, $p->avg_length, $p->max_length, $p->standard_deviations], [7156780, 364995780, '51.00', 51, ['0.00', 0.02]], 'header of fastqcheck file was parsed correctly';
$p->next_record;
my $first_record = [@{$p->parsed_record}];
my ($bases, $quals) = $p->avg_base_quals();
my $num_records = 1;
while ($p->next_record) {
    $num_records++;
}
is_deeply [$num_records, $first_record, $quals, $p->avg_qual],
                                        [52,
                                         [qw(0 25.4 22.0 20.4 28.7  3.6    0   0   0   0   0  36  19  26
                                             3 1  16  16   5   7  19  29  10  13  17  17  21  18  17  20  16
                                             23  22  20  18  23 20  23  27  27  31  41  91 114 115  40  20
                                             15.2)],
                                         [qw(33.6929292929293 32.4489383215369 33.8459214501511 32.920282542886 33.148743718593 33.204843592331 33.832995951417 34.5202429149798 33.6373737373737
                                             33.9384460141271 34.5242914979757 33.6814964610718 34.1830131445905 35.34375 34.0594758064516 31.3387259858443 32.8558467741936 32.8528225806452
                                             29.8375378405651 32.7464646464646 32.8977732793522 32.4495967741936 32.4118831822759 32.8042381432896 33.0604838709677 32.7993951612903
                                             31.6368209255533 32.3256048387097 22.9908814589666 22.5510616784631 24.2992922143579 22.5556680161943 21.9604863221885 25.3313131313131
                                             24.2096774193548 25.0443101711984 24.1553985872856 23.6963562753036 23.314459049545 21.7979797979798 22.3222222222222 19.2127016129032
                                             20.8827098078868 20.7700101317123 19.837044534413 18.9494438827098 18.0485829959514 17.7537993920973 16.8546922300706 17.1151515151515 17.3699596774194)],
                                         27.8919469928644], 'body of fastqcheck file was parsed correctly';

# fasta
$p = VRPipe::Parser->create('fasta', {file => file(qw(t data S_suis_P17.fa))});
@records = ();
my $pr = $p->parsed_record;
while ($p->next_record) {
    push(@records, $pr->[0].'+'.length($pr->[1]));
}
is_deeply [@records], ["fake_chr1+290640", "fake_chr2+1716851"], 'fasta file was parsed correctly';

# dict
$p = VRPipe::Parser->create('dict', {file => file(qw(t data S_suis_P17.fa.dict))});
@records = ();
$pr = $p->parsed_record;
while ($p->next_record) {
    push(@records, {%{$pr}});
}
my @common_dict = (UR => 'ftp://s.suis.com/ref.fa', AS => 'SSuis1', SP => 'S.Suis');
is_deeply [@records], [{SN => 'fake_chr1', LN => 290640, M5 => '55f9584cf1f4194f13bbdc0167e0a05f', @common_dict},
                       {SN => 'fake_chr2', LN => 1716851, M5 => '6dd2836053e5c4bd14ad49b5b2f2eb88', @common_dict}], 'dict file was parsed correctly';
is $p->total_length, 2007491, 'total_length for dict files worked';

# bam (we have Parser_bam.t for more in-depth testing of bam parsing)
$p = VRPipe::Parser->create('bam', {file => file(qw(t data file.bam))});
is $p->sam_version, '1.0', 'sam version could be parsed from bam file';
is $p->sequence_info('2', 'LN'), 243199373, 'chrom length could be parsed from bam file';
is $p->program_info('bwa', 'VN'), '0.5.5', 'program info could be parsed from bam file';
is $p->readgroup_info('SRR035022', 'LB'), 'Solexa-16652', 'readgroup info could be parsed from bam file';
is_deeply [$p->samples], ['NA06984'], 'could get all samples from bam file';
$num_records = 0;
while ($p->next_record) {
    $num_records++;
}
is $num_records, 3, 'correct number of records found in bam file';

# bamcheck
$p = VRPipe::Parser->create('bamcheck', {file => file(qw(t data parser.bamcheck))});
is $p->sequences, 2000, 'bamcheck sequences method worked';
is $p->is_paired, 1, 'bamcheck is_paired worked';
is $p->bases_mapped_cigar, 58663, 'bamcheck bases_mapped_cigar worked';
is $p->error_rate, '2.151271e-02', 'bamcheck error_rate worked';
is $p->insert_size_standard_deviation, 70.9, 'bamcheck insert_size_standard_deviation worked';
# some (newer?) bamcheck files have more lines than others:
$p = VRPipe::Parser->create('bamcheck', {file => file(qw(t data parser.bamcheck_mq0))});
is $p->reads_mq0, 1859041, 'reads_mq0 worked';
is $p->pairs_with_other_orientation, 5857, 'pairs_with_other_orientation worked';
# use a bamcheck with a simple COV section to test coverage methods
$p = VRPipe::Parser->create('bamcheck', {file => file(qw(t data parser.bamcheck_simplecov))});
is_deeply [$p->coverage(0), $p->coverage(1), $p->coverage(2), $p->coverage(3), $p->coverage(4), $p->coverage(5), $p->coverage(6), $p->coverage(7)], [0, 27893, 6485, 1188, 355, 49, 49, 0], 'coverage worked';
is_deeply [$p->cumulative_coverage(1), $p->cumulative_coverage(2), $p->cumulative_coverage(3), $p->cumulative_coverage(4), $p->cumulative_coverage(5), $p->cumulative_coverage(6), $p->cumulative_coverage(7)], [36019, 8126, 1641, 453, 98, 49, 0], 'cumulative_coverage worked';
is $p->mean_coverage, 1.29, 'mean_coverage worked';
# and check that some of the other section methods work
is_deeply $p->gc_depth(), [[qw(0.0	28.571	0.000	0.000	0.000	0.000	0.000)],
                           [qw(0.4	42.857	0.238	0.238	0.238	0.238	0.238)],
                           [qw(37.0	71.429	0.381	0.381	0.381	0.383	0.383)],
                           [qw(38.0	85.714	0.402	0.402	0.402	0.402	0.402)],
                           [qw(39.0	100.000	0.419	0.419	0.419	0.419	0.419)]], 'gc_depth worked';
is_deeply $p->read_lengths(), [[54, 862]], 'read_lengths worked';
is_deeply $p->indel_dist(), [[1, 5, 6], [3, 0, 1]], 'indel_dist worked';

# bas
$p = VRPipe::Parser->create('bas', {file => file(qw(t data example2.bas))});
$pr = $p->parsed_record;
@records = ();
while ($p->next_record) {
    push(@records, $pr->[0].'=='.$pr->[20]);
}
is_deeply \@records, ['NA00001.ILLUMINA.bwa.unknown_population.unknown_analysisgroup.20100208==122',
                      'NA00002.ILLUMINA.bwa.unknown_population.unknown_analysisgroup.20100208==122',
                      'NA00003.ILLUMINA.bwa.unknown_population.unknown_analysisgroup.20100208==122'], 'correct number of records found in bas file, and seems to be fully parsed';
is $p->total_bases, 345000, 'total_bases worked';
is $p->duplicate_bases, 366, 'duplicate_bases worked';

# fastq
{
    my $fq_file = Path::Class::File->new(qw(t data parser.fastq));
    my $gz_file =  Path::Class::File->new(qw(t data parser.fastq.gz));
    
    my $fqp = VRPipe::Parser->create('fastq', {file => $fq_file});
    my $pr = $fqp->parsed_record;
    
    # parse the first sequence
    $fqp->next_record;
    is_deeply $pr, ['SRR001629.5', 'GGGGGCAATGCTGGGGTCGTCCTCCTCAACTCGCTCCAGGGGCCAGGGGATACCGCTCATATCACTAAGGGCGGTGCCCAGGTAGAGGAGCTCGCGATAGTCCCATTCAATGGACGTGTACCGGATGTTTAGGAGAGGCAGGGAGGCGATGATCTGGCATGTGTGCCGCAGGTGTGTCAGGAGGTCGTCAA', '88888>>BBBB>A@@@@BBBBBBBBBAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB@@@BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB@@@ABBBBBBBB'], 'parsed data for first sequence';
    
    # parse the last line
    while ($fqp->next_record) { next; };
    is_deeply $pr, ['SRR001629.2599', 'CACGTGTCCTCAACCTTGGCAAAATAAACTTTCAAAATTAACTGAGACCTATCTCAGATTTTCGGGGTTCACAGTAGCAAGGAAGTGGGGTCTGAGAGACGCCCC', 'BBBBBBBBBBBBBBBBBBAAAAAA?BBB?AAA?>>>>@@@@?B??BBBBBBBBBBBBB@@@@AAB@@@@AABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB'], 'parsed data for last sequence';
    
    # needs to work on .gz files as well
    $fqp = VRPipe::Parser->create('fastq', {file => $gz_file});
    $pr = $fqp->parsed_record;
    while ($fqp->next_record) { next; };
    is $pr->[1], 'CACGTGTCCTCAACCTTGGCAAAATAAACTTTCAAAATTAACTGAGACCTATCTCAGATTTTCGGGGTTCACAGTAGCAAGGAAGTGGGGTCTGAGAGACGCCCC', 'seq correct for last read of gz compressed fastq ';
    
    # qual to ints
    is_deeply [$fqp->qual_to_ints('!"#$%&\'()*+,5?DIS]~')], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 20, 30, 35, 40, 50, 60, 93], 'qual_to_ints test worked';
}

# sequence.index
{
    my $sip = VRPipe::Parser->create('sequence_index', {file => file(qw(t data parser.sequence_index))});
    my $rh = $sip->parsed_record;
    
    my @expected_data = (['data/NA19238/sequence_read/ERR000018.recal.fastq.gz',
                          'ee05c1a260621d8840ddf3028ebb2355',
                          'ERR000018',
                          'SRP000032',
                          '1000Genomes Project Pilot 2',
                          'BGI',
                          'ERA000013',
                          '',
                          'SRS000212',
                          'NA19238',
                          'YRI_1',
                          'ERX000014',
                          'ILLUMINA',
                          'Illumina Genome Analyzer',
                          'HU1000RADCAASE',
                          'BGI-FC307N0AAXX',
                          'BGI-FC307N0AAXX_5',
                          '',
                          'SINGLE',
                          '',
                          '0',
                          '',
                          '',
                          '9612363',
                          '346045068'],
                         ['data/NA12282/sequence_read/SRR015438_2.recal.fastq.gz',
                          '80d7ee75e062bbd76756df1dc94c6539',
                          'SRR015438',
                          'SRP000033',
                          '1000Genomes Project Pilot 3',
                          'WUGSC',
                          'SRA008537',
                          '',
                          'SRS000619',
                          'NA12282',
                          'CEPH - 2',
                          'SRX004024',
                          'ILLUMINA',
                          'Illumina Genome Analyzer II',
                          '2773138721',
                          '28075',
                          'HWI-EAS289_3150M',
                          '260',
                          'PAIRED',
                          'data/NA12282/sequence_read/SRR015438_1.recal.fastq.gz',
                          '0',
                          '',
                          '',
                          '12938797',
                          '659878647']);
    
    # parse the first line
    $sip->next_record;
    my $expected = shift @expected_data;
    is_deeply $rh, $expected, 'parsed data for first line of a sequence.index';
    
    # get info on a particular lane from line 10
    is $sip->lane_info('ERR000025', 'sample_name'), 'NA19240', 'lane_info test when not yet reached';
    is $rh->[0], $expected->[0], 'using lane_info doesn\'t change our result holder';
    $sip->next_record;
    is $rh->[0], 'data/NA19238/sequence_read/ERR000019.recal.fastq.gz', 'using lane_info doesn\'t mess with next_result';
    
    # get info on a particular lane from line 5
    is $sip->lane_info('ERR000020', 'fastq_file'), 'data/NA19240/sequence_read/ERR000020_2.recal.fastq.gz', 'lane_info test when allready seen';
    
    # parse the last line
    while ($sip->next_record) { next; };
    $expected = shift @expected_data;
    is_deeply $rh, $expected, 'parsed data for last line';
    
    # test get_lanes() on headed
    my @all_lanes = $sip->get_lanes;
    is $all_lanes[0], 'ERR000018', 'got first lane with get_lanes on headed file';
    my %all_lanes = map { $_ => 1 } @all_lanes;
    ok ! defined $all_lanes{RUN_ID}, 'RUN_ID did not get treated as a lane';
    
    # try parsing a sequence.index with no header line
    $sip = VRPipe::Parser->create('sequence_index', {file => file(qw(t data parser.sequence_index_headerless))});
    $rh = $sip->parsed_record;
    
    is $sip->lane_info('ERR000044', 'sample_name'), 'NA18550', 'headerless file parse worked';
    $sip->next_record;
    is $rh->[0], 'data/NA18550/sequence_read/ERR000044_1.recal.fastq.gz', 'got first line correctly';
    
    # test get_lanes() on headerless
    @all_lanes = $sip->get_lanes;
    is $all_lanes[0], 'ERR000044', 'got first lane with get_lanes on headerless file';
    is $all_lanes[-1], 'SRR014220', 'got last lane';
    is @all_lanes, 8723, 'got all lanes';
    is $rh->[0], 'data/NA18550/sequence_read/ERR000044_1.recal.fastq.gz', 'getting all lanes didn\'t alter our result holder';
    while ($sip->next_record) {
        next;
    }
    is $rh->[0], 'data/NA19093/sequence_read/SRR014220_2.recal.fastq.gz', 'got last line correctly';
    
    my @lanes = $sip->get_lanes(sample_name => 'NA11994');
    is $lanes[0], 'SRR003428', 'got first lane with a given sample_name';
    is $lanes[-1], 'SRR014158', 'got last lane with a given sample_name';
    is @lanes, 120, 'got all lanes with a given sample_name';
    
    @lanes = $sip->get_lanes(ignore_withdrawn => 1);
    is @lanes, 8508, 'got all non-withdrawn lanes';
    @lanes = $sip->get_lanes(ignore_INSTRUMENT_PLATFORM => 'solid');
    is @lanes, 4544, 'got all non-solid lanes';
    
    # problem run_id where it is both withdrawn and not withdrawn
    my @answers = $sip->lane_info('ERR000061', 'withdrawn');
    is_deeply \@answers, [0, 1], 'knew that a given lane was both withdrawn and not';
    is $sip->lane_info('ERR000061', 'withdrawn'), 1, 'knew that the lane was noted as withdrawn more times than not';
    
    # in 2010, sequence.index format changed by adding a new ANALYSIS_GROUP column
    $sip = VRPipe::Parser->create('sequence_index', {file => file(qw(t data parser.sequence_index_2010))});
    is $sip->lane_info('ERR000018', 'analysis_group'), 'high coverage', 'ANALYSIS_GROUP is parsable';
}

exit;