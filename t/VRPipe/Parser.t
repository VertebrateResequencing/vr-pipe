#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class qw(file);

BEGIN {
    use Test::Most tests => 16;
    
    use_ok('VRPipe::Parser');
    
    use TestPersistentReal;
}

ok my $p = VRPipe::Parser->create('lsf', {file => file(qw(t data lsf.stdout))}), 'could create an lsf parser';
$p->fh;
is_deeply [$p->file, $p->_vrpipe_file->path], [file(qw(t data lsf.stdout)), file(qw(t data lsf.stdout))->absolute], 'file details are correct';

throws_ok {$p = VRPipe::Parser->create('foo', {});} qr/Invalid implementation class/, 'throws when asked to create an invalid parser';

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

throws_ok {$p = VRPipe::Parser->create('cat', {file => 't/data/not.extant'}); $p->next_record; } qr/does not exist, so cannot be opened for reading/, 'throws when asked to parse a non-existant file';

undef $p;
$p = VRPipe::Parser->create('cat', {file => file(qw(t data file.cat))});
my @records;
while ($p->next_record) {
    push(@records, join("\n", @{$p->parsed_record}));
}
is_deeply [@records], ["first line of 4th record\nsecond line of 4th record",
                       "",
                       "first line of 2nd record\nsecond line of 2nd record",
                       "first line of 1st record\nsecond line of 1st record"], 'cat file was parsed correctly';

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

exit;