#!/usr/bin/env perl
use strict;
use warnings;
use File::Temp qw/tempdir/;
use File::Spec;
use File::Path qw/make_path/;
use Class::Unload;

BEGIN {
    #use Test::Expect;
    use Test::Most tests => 53;
    
    use_ok('VRPipe::Config');
    use_ok('t::VRPipe::Config');
    use_ok('VRPipe::Base::Configuration::Env');
}

my $tmp_dir = tempdir(CLEANUP => 1);
my $temp_path = File::Spec->catdir($tmp_dir, qw(VRPipe));
make_path($temp_path);

unshift(@INC, $tmp_dir);
ok my $vrp_config = reload(), 'made a new Config object';
is $vrp_config->config_module_path, File::Spec->catfile($temp_path, 'TestConfig.pm'), 'config_module_path is as expected';

ok my $option = $vrp_config->next_option, 'got an option from next_option';
isa_ok $option, 'VRPipe::Base::Configuration::Option';
is $option->question,        'First question?', 'option has correct question';
is $option->question_number, 1,                 'also has a question number';
is $option->key,             'first_key',       'key is correct';
is $option->value, $vrp_config->first_key, '$option->key is the same as $config->$key_name()';
is $option->value, 'foo', 'value defaults to its specified default';
is $option->value('bar'),   'bar', 'value could be set to bar';
is $option->value('undef'), undef, 'value can be set to string undef to return actual undef';
is_deeply $option->valid, [qw(foo bar make_second_skip)], 'valid contains allowed values';
throws_ok { $option->value('baz') } qr/not a valid value/, 'value throws if supplied an invalid value';

ok $option = $vrp_config->next_option, 'got an option from next_option';
is $option->question, 'Second question?', 'second option has correct question';
is $option->value, undef, 'value defaults to undef';
{
    local $ENV{vrpipe_second_key} = 'foo';
    is $option->value, undef, 'value default does not change with env vars after first accessed';
    $vrp_config = reload();
    ok $option = $vrp_config->option(number => 2), 'option could be retreived by number';
    is $option->value, 'foo', 'value defaults to an env var if defined prior to first access';
    is $option->value('bar'), 'bar', 'value can still be set directly';
    ok $vrp_config->write_config_module, 'able to write_config_module';
    $vrp_config = reload();
    is $vrp_config->second_key, 'bar', 'second_key set value was retained after reloading';
    is $vrp_config->second_key(undef), undef, 'second_key can be cleared and becomes undef';
    $vrp_config->write_config_module;
    $vrp_config = reload();
    is $vrp_config->second_key, 'foo', 'after reload, second_key defaults back to env var having been undefined';
}

$vrp_config = reload();
$vrp_config->next_option;
$vrp_config->next_option;
ok $option = $vrp_config->next_option, 'got a third option from next_option';
is $option->question, 'third_key?', 'option has correct question';
is $option->value,    undef,        'value defaults to undef';
{
    local $ENV{TVRPIPE_THIRDKEY} = 'foo';
    $vrp_config = reload();
    ok $option = $vrp_config->option(key => 'third_key'), 'option could be retreived by key';
    is $option->value, 'foo', 'value defaults to a config-defined special env var';
    local $ENV{vrpipe_test_var} = 'bar';
    my $env_obj = VRPipe::Base::Configuration::Env->new(variable => 'vrpipe_test_var');
    my $return_val = $vrp_config->third_key($env_obj);
    is_deeply [ref($return_val), "$return_val"], ['VRPipe::Base::Configuration::Env', 'bar'], 'third_key could be set with an Env object and returns the env var';
    local $ENV{vrpipe_test_var2} = 'baz';
    ok $option->env('vrpipe_test_var2'), 'could set $option->env';
    $return_val = $option->value;
    is_deeply [ref($return_val), "$return_val"], ['VRPipe::Base::Configuration::Env', 'baz'], 'value could be set via env and returns the env var';
    $vrp_config->write_config_module;
    $vrp_config = reload();
    is $vrp_config->third_key, 'baz', 'third_key returns the custom env var value after reload';
    delete $ENV{vrpipe_test_var2};
    $vrp_config = reload();
    $return_val = $vrp_config->third_key;
    is_deeply [$return_val->variable, "$return_val"], ['vrpipe_test_var2', ''], 'when the env var is undefined, we still return the Env object, but the value is empty string';
    $vrp_config->third_key(undef);
    $vrp_config->write_config_module;
}

# skipping questions and changing defaults based on past answers
$vrp_config = reload();
$option     = $vrp_config->next_option;
is $option->key, 'first_key', 'after a reload we start with the first option again...';
$option->value('make_second_skip');
$option = $vrp_config->next_option;
is $option->key, 'third_key', '... but the second option gets skipped based on value of first option';
is $option->value, 'based on skipped second key', 'the third option got a special default based on second option being skipped';

# prompts
# *** tried using Test::Expect to test prompt(), but couldn't get it to work
#     beyond the first test.
#expect_run(
#    command => "perl t/VRPipe/Config.pl",
#    prompt  => 'First question? <foo|bar|make_second_skip> [foo] ',
#    quit    => 'quit',
#);
#expect 'ping', "'ping' was not a valid answer for that question; try again:", 'prompt supplied invalid answer prompts again';
#expect 'foo', 'Second question? [] ', 'second question gets asked after first question answered correctly';
#expect 'bar', 'third_key? [] ', 'third question gets asked after second answer';
#expect_send 'baz', 'exits after third answer';

# security
$vrp_config = reload(1);
my $injected_first_option = $vrp_config->next_option;
is $injected_first_option->question_number, 0,                                                                                                                                   'when we have secure options, a 0th option is injected';
is $injected_first_option->key,             'encryption_key_file',                                                                                                               'the 0th option has the correct key';
is $injected_first_option->question,        'Passwords you enter will be encrypted; where should your encryption key be stored? (it is up to you to properly secure this file)', 'the 0th option asks about the key file';
is $injected_first_option->value,           undef,                                                                                                                               'key file has no default';
my $secret_option = $vrp_config->option(key => 'secret_key');
ok $secret_option->secure, 'got the secret_key option which advertises itself as secure';
is $secret_option->value('foo'), 'foo', 'a secret option can be set without a key file...';
throws_ok { $vrp_config->write_config_module() } qr/Secure config options cannot be set without the encryption_key_file option being set/, "... but you can't write a secure value when the key file has not been defined";
my $keyfile = File::Spec->catfile($tmp_dir, 'keyfile');
is $injected_first_option->value($keyfile), $keyfile, 'the key file location could be set';
ok $vrp_config->write_config_module, 'once the key file has been defined, we can write a secure value';
$vrp_config = reload();
eval "require VRPipe::TestConfig;";
my $raw_config_hash = VRPipe::TestConfig->get_config;
like $raw_config_hash->{secret_key}, qr/^Salted__/, 'the secure value is stored encrypted';
$secret_option = $vrp_config->option(key => 'secret_key');
is $secret_option->value, 'foo', 'yet the secure value can still be retrieved unencrypted after a reload';
$secret_option->env('vrpipe_securetest_env');
$vrp_config->write_config_module;
$vrp_config = reload();
my $securetest_env = $vrp_config->secret_key;
isa_ok $securetest_env, 'VRPipe::Base::Configuration::Env', 'an env assigned to a secure value comes back an env, not an encrypted mess;';

# just check that looping with next_option works correctly
$vrp_config = reload(1);
my $options = 0;
while (my $test_option = $vrp_config->next_option) {
    $options++;
}
is $options, 5, 'got correct number of options';

exit;

sub reload {
    Class::Unload->unload('VRPipe::TestConfig');
    my $c = t::VRPipe::Config->new(config_module => 'VRPipe::TestConfig');
    $c->next_option unless shift;
    return $c;
}
