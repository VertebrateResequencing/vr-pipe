#!/usr/bin/env perl
use strict;
use warnings;
use File::Temp qw/tempdir/;
use File::Spec;
use File::Path qw/make_path/;
use Class::Unload;

BEGIN {
    use Test::Most tests => 29;
    
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
my @options = $vrp_config->get_options;
is @options, 3, 'got correct number of options';

my $option = $options[0];
is ref($option), 'HASH', 'an option is a hashref';
is $option->{question}, 'First question?', 'option has correct question';
is $option->{question_number}, 1, 'also has a question number';
is $option->{key}, 'first_key', 'key is correct';
is $vrp_config->first_key, 'foo', 'first_key defaults to its specified default';
is $vrp_config->first_key('bar'), 'bar', 'first_key could be set to bar';
is_deeply $option->{valid}, [qw(foo bar)], 'valid contains allowed values';
is $vrp_config->first_key('baz'), 'baz', 'first_key can be set to an invalid value';

$option = $options[1];
is $option->{question}, 'Second question?', 'second option has correct question';
is $vrp_config->second_key, undef, 'second_key defaults to undef';
{
    local $ENV{vrpipe_second_key} = 'foo';
    is $vrp_config->second_key, undef, 'second_key default does not change with env vars after first accessed';
    $vrp_config = reload();
    is $vrp_config->second_key, 'foo', 'second_key defaults to an env var if defined prior to first access';
    is $vrp_config->second_key('bar'), 'bar', 'second_key can still be set directly';
    ok $vrp_config->write_config_module, 'able to write_config_module';
    $vrp_config = reload();
    is $vrp_config->second_key, 'bar', 'second_key set value was retained after reloading';
    is $vrp_config->second_key(undef), undef, 'second_key can be cleared and becomes undef';
    $vrp_config->write_config_module;
    $vrp_config = reload();
    is $vrp_config->second_key, 'foo', 'after reload, second_key defaults back to env var having been undefined';
}

$option = $options[2];
is $option->{question}, 'Third question?', 'option has correct question';
is $vrp_config->third_key, undef, 'third_key defaults to undef';
{
    local $ENV{TVRPIPE_THIRDKEY} = 'foo';
    $vrp_config = reload();
    is $vrp_config->third_key, 'foo', 'third_key defaults to a config-defined special env var';
    local $ENV{vrpipe_test_var} = 'bar';
    my $env_obj = VRPipe::Base::Configuration::Env->new(variable => 'vrpipe_test_var');
    my $return_val = $vrp_config->third_key($env_obj);
    is_deeply [ref($return_val), "$return_val"], ['VRPipe::Base::Configuration::Env', 'bar'], 'third_key could be set with an Env object and returns the env var';
    $vrp_config->write_config_module;
    $vrp_config = reload();
    is $vrp_config->third_key, 'bar', 'third_key returns the custom env var value after reload';
    delete $ENV{vrpipe_test_var};
    $vrp_config = reload();
    $return_val = $vrp_config->third_key($env_obj);
    is_deeply [ref($return_val), "$return_val"], ['VRPipe::Base::Configuration::Env', ''], 'when the env var is undefined, we still return the Env object, but the value is empty string';
}

exit;

sub reload {
    Class::Unload->unload('VRPipe::TestConfig');
    return t::VRPipe::Config->new(config_module => 'VRPipe::TestConfig');
}