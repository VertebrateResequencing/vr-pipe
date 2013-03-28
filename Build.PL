#!/usr/bin/env perl
# Build script mainly used to run tests.
#
# Do this once, first:
# perl Build.PL
# Then before each commit, run the test for the module you altered:
# ./Build test --verbose --test_files t/appropriatetestscript.t
# And if you did something major, run the whole test suite:
# ./Build test
#
# less often, when interested:
# ./Build testcover
# ./Build testpod
# ./Build testpodcoverage
#
# setenv DEVEL_COVER_OPTIONS '+ignore,5.8.8'
# beforehand if coverage includes extraneous stuff
#

use strict;
use warnings;
use lib 'build_lib';
use VRPipeBuild;

my $build = VRPipeBuild->new(
    module_name          => 'VRPipe',
    dist_version         => 0.148,
    dist_author          => 'Vertebrate Resequencing group at the Sanger Institute',
    dist_abstract        => 'A system for setting up, running and tracking many jobs in pipelines on a compute farm.',
    license              => 'gpl',
    recursive_test_files => 1,
    pm_files             => get_pm_files,
    script_files         => 'scripts',
    requires             => required_modules
);

$build->create_site_config;

# *** we need _build/build_params to contain the SiteConfig.pm file we just
# created in our modules directory, but I can't figure out how to update that
# file. Without the following hack, our SiteConfig.pm won't get installed with
# ./Build install unless 'perl Build.PL' is run once again
$build->{properties}->{pm_files} = get_pm_files;

$build->create_build_script;

exit;
