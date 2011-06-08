#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class qw(file);

BEGIN {
    use Test::Most tests => 12;
    
    use_ok('VRPipe::Parser');
    
    use TestPersistentReal;
}

ok my $p = VRPipe::Parser->create('lsf', {file => 't/data/lsf.stdout'}), 'could create an lsf parser';
$p->fh;
is_deeply [$p->file, $p->_vrpipe_file->path], ['t/data/lsf.stdout', file('t/data/lsf.stdout')->absolute], 'file details are correct';

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
                              'normal'], 'parsed_record contains all the correct details';

ok $p->next_record, 'next_record worked again';
is $p->cpu_time, 5.05, 'cpu_time worked after next_record calls';

ok ! $p->next_record, 'next_record returns false when no more records';

undef $p;
$p = VRPipe::Parser->create('lsf', {file => 't/data/lsf.stdout'});
$p->next_record;
is $p->parsed_record->[4], 4.75, 'next_record and parsed_record work on the first (last) record';


exit;