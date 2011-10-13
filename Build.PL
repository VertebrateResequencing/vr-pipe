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
    dist_version         => 0.12,
    dist_author          => 'Vertebrate Resequencing group at the Sanger Institute',
    dist_abstract        => 'A system for setting up, running and tracking many jobs in pipelines on a compute farm.',
    license              => 'perl',
    recursive_test_files => 1,
    pm_files             => get_pm_files(),
    script_files         => 'scripts',
    requires             => { perl => '5.8.8',
                              'B::Deparse' => 0,
                              'Crypt::CBC' => 0,
                              'Cwd' => 0,
                              'Data::Compare' => 0,
                              'Data::Dumper' => 0,
                              'DateTime' => 0,
                              'DBIx::Class' => 0,
                              'DBIx::Class::DeploymentHandler' => 0,
                              'Devel::GlobalDestruction' => 0,
                              'Digest::MD5' => 0,
                              'File::Copy' => 0,
                              'File::Fetch' => 0,
                              'File::HomeDir' => 0,
                              'File::ReadBackwards' => 0,
                              'File::Spec' => 0,
                              'File::Temp' => 0,
                              'IO::Capture::Stderr' => 0,
                              'IO::Uncompress::AnyUncompress' => 0,
                              'Moose' => 0,
                              'MooseX::AbstractFactory' => 0,
                              'MooseX::Aliases' => 0,
                              'MooseX::Declare' => 0,
                              'Net::FTP::Robust' => 0,
                              'Net::SSH' => 0,
                              'Parallel::ForkManager' => 0,
                              'Path::Class' => 0,
                              'POSIX' => 0,
                              'Sys::CPU' => 0,
                              'Sys::Hostname' => 0,
                              'Time::Format' => 0 } );

$build->create_site_config;
$build->create_build_script;

exit;

sub get_pm_files {
	my %pm_files;
	foreach my $module (check_dir('modules')) {
		my $in_lib = $module;
		$in_lib =~ s/^modules/lib/;
		$pm_files{$module} = $in_lib;
	}
	return \%pm_files;
}

sub check_dir {
        my $dir = shift;
        opendir(my $dir_handle, $dir);
	
	my @pm_files;
        foreach my $thing (readdir($dir_handle)) {
                if ($thing =~ /^\.+$/) { next; }
                $thing = $dir."/".$thing;

                if (-d $thing) {
                        push(@pm_files, check_dir($thing));
                        next;
                }
		
		if ($thing =~ /\.pm$/) {
			push(@pm_files, $thing);
		}
        }
	return @pm_files;
}
