#!/usr/bin/env perl
use strict;
use warnings;

use File::Temp qw/tempdir/;
use File::Path qw/make_path/;

use t::VRPipe::Config;

my $tmp_dir = tempdir(CLEANUP => 1);
my $temp_path = File::Spec->catdir($tmp_dir, qw(VRPipe));
make_path($temp_path);

unshift(@INC, $tmp_dir);
my $vrp_config = t::VRPipe::Config->new(config_module => 'VRPipe::TestConfig');

while (my $option = $vrp_config->next_option) {
    $option->prompt;
}

$vrp_config->write_config_module();

exit;
