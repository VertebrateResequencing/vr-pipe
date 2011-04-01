#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 49;
    
    use_ok('VRPipe::Base');
    use_ok('File::Spec');
}

package VRPipe::Test;
use VRPipe::Base;
class VRPipe::Test {
    has 'foo' => (
        is  => 'rw',
        isa => 'Str',
    );
    
    has 'file' => (
        is  => 'rw',
        isa => File,
        coerce => 1
    );
    
    has 'optional_file' => (
        is  => 'rw',
        isa => MaybeFile,
        coerce => 1
    );
}
1;

package VRPipe::Test2;
use VRPipe::Base;
class VRPipe::Test2 {
    has 'foogu' => (
        is  => 'rw',
        isa => 'Str',
    );
}
1;

package main;
my $base = VRPipe::Test->new(foo => 'bar');
isa_ok $base, 'VRPipe::Test';
is $base->foo, 'bar', 'args sent to new are set in the corresponding method';

# verbose and warn testing
is $base->verbose, 0, 'default verbose value';
warning_is { $base->warn('simple msg') } { carped => 'simple msg' }, 'verbose 0 warning is simple carp';
is $base->verbose(-1), -1, 'verbose -1 could be set';
warning_is { $base->warn('simple msg') } '', 'no warning when verbose -1';
is $base->verbose(1), 1, 'verbose 1 could be set';
warning_is { $base->warn('simple msg') } { carped => 'simple msg' }, 'verbose 1 warning is a cluck';
is $base->verbose(2), 2, 'verbose 2 could be set';
throws_ok { $base->warn('thrown msg') } qr/thrown msg/, 'verbose 2 warning is a throw';
is $base->verbose(0), 0, 'verbose 0 could be set';

$base->verbose(1);
my $base2 = VRPipe::Test->new(foo => 'bar', verbose => 2);
is $base2->verbose, 2, 'verbose can be set via new';
is $base->verbose, 1, 'verbose on an instance is not set for the whole class';
my $v = VRPipe::Test->verbose(-1);
is $v, -1, 'verbose could be set on the class';
is $base->verbose, -1, 'the class set affects instances';
my $other = VRPipe::Test2->new(foogu => 'bazgu');
is $other->verbose, -1, 'the class set affects instances of other classes';
is $base->verbose(2), -1, 'instance set does not override global';
VRPipe::Test2->clear_verbose();
ok 1, 'global verbose can be unset';
is $other->verbose(1), 1, 'instance set works again after unsetting global';
is $base->verbose(), 2, 'instance sets whilst global override in place were not forgotten';

# throw and debug
throws_ok { $base->throw('thrown msg') } qr/thrown msg/, 'throw works';
warning_like { $base->debug('debug msg') } qr/debug msg/, 'debug message when verbose > 0';
$base->verbose(0);
warning_is { $base->debug('debug msg') } '', 'no debug message when verbose <= 0';

# type constraints
my $fname = 'sdf';
is $base->file($fname), $fname, "'$fname' passed the file constraint";
$fname = '|!"£$%^&*()_+{}~@:<>?';
throws_ok { $base->file($fname) } qr/does not seem like a file/, "'$fname' fails the file constraint";
$fname = 'foo/bar/baz.txt';
my $fobj = $base->file($fname)->as_foreign('Win32');
is "$fobj", 'foo\bar\baz.txt', "'$fname' passed the file constraint (converting to $fobj for Win32)";
is $fobj->basename, 'baz.txt', 'file() returned an object we could get the correct basename from';
$fname = [qw(foo bar baz.txt)];
$fobj = $base->file($fname);
like "$fobj", qr/^foo.bar.baz\.txt$/, "[@{$fname}] as array ref passed the file constraint (converting to $fobj for this platform)";
is $fobj->basename, 'baz.txt', 'basename also works given array ref input';
throws_ok { $base->file(undef) } qr/no file specified/, "undef fails the file constraint";
is $base->optional_file(undef), undef, 'undef passes the MaybeFile constraint';

# register_for_cleanup difficult to test for properly...
can_ok $base, qw(register_for_cleanup);
$base->register_for_cleanup('foo');
$base->register_for_cleanup('verbose');
$base->register_for_cleanup('not_correct');
is_deeply {$base->cleanup_methods}, {foo => 1, verbose => 1}, 'existing methods can be registered for cleanup';
$base->unregister_for_cleanup('verbose');
is_deeply {$base->cleanup_methods}, {foo => 1}, 'method can be unregistered for cleanup';

# unlinking can be tested, alongside logging
my $tfile1 = File::Spec->catfile('t', 'data', 'VRPipe_base_unlink_test_1');
my $tfile2 = File::Spec->catfile('t', 'data', 'VRPipe_base_unlink_test_2');
system("touch $tfile1; touch $tfile2");
$base->register_for_unlinking($tfile1, $tfile2);
is_deeply {$base->files_to_unlink}, {$tfile1 => 1, $tfile2 => 1}, 'files can be registered for unlinking';
$base->unregister_for_unlinking($tfile2);
is_deeply {$base->files_to_unlink}, {$tfile1 => 1}, 'files can be unregistered for unlinking';
$base->register_for_unlinking($tfile2);
ok -e $tfile1, 'file1 exists prior to unlink attept';
ok -e $tfile2, 'file2 exists prior to unlink attept';
my $default_log_file = $base->log_file();
like $default_log_file, qr/VRPipe.log$/, 'default log file has correct name';
is $base->log_file($tfile1), $tfile1, 'log file location could be changed';
is $base->write_logs, 0, 'logging off by default';
$base->log('message');
ok ! -s $tfile1, 'logging while logging is off does nothing';
is $base->write_logs(1), 1, 'logging could be turned on';
$base->log('message');
my $log_size = -s $tfile1;
cmp_ok $log_size, '>', 0, 'logging while logging is on does something';
$base->log('message2');
my $log_size2 = -s $tfile1;
cmp_ok $log_size2, '>', $log_size, 'adding a log message appends';
undef $base;
ok ! -e $tfile1, 'file1 unlinked after object destruction';
ok ! -e $tfile2, 'file1 unlinked after object destruction';

exit;
