
=head1 NAME

VRPipeBuild - Module::Build subclass with VRPipe-specific bits

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

package VRPipeBuild;

BEGIN { unshift(@INC, './modules') }
use base qw(Module::Build Exporter);
@EXPORT = qw(get_pm_files required_modules);
use Config;
use File::Spec;

# we're in the odd position of wanting to call the code we're trying to install
# during the install process. If we're missing module dependencies
# create_site_config() below will fail because it calls VRPipe::Config and then
# the Build script won't get created, which means the user can't use it to
# auto-install our dependencies! So create_site_config needs to check we have
# all deps installed.

sub required_modules {
    return {
        perl                             => '5.8.8',
        'AnyEvent'                       => 0,
        'AnyEvent::ForkManager'          => 0,
        'B::Deparse'                     => 0,
        'Class::Unload'                  => 0,
        'Continuity'                     => 0,
        'Crypt::CBC'                     => 0,
        'Crypt::Blowfish'                => 0,
        'Crypt::Random'                  => 0,
        'Cwd'                            => 0,
        'Data::Compare'                  => 0,
        'Data::Dumper'                   => 0,
        'DateTime'                       => 0,
        'DateTime::Format::Natural'      => 0,
        'DBI'                            => 0,
        'DBIx::Class'                    => '>= 0.08200',
        'DBIx::Class::DeploymentHandler' => 0,
        'Devel::GlobalDestruction'       => 0,
        'Digest::MD5'                    => 0,
        'Email::Sender::Simple'          => 0,
        'Email::Simple::Creator'         => 0,
        'EV'                             => 0,
        'File::Copy'                     => 0,
        'File::Fetch'                    => 0,
        'File::HomeDir'                  => 0,
        'File::Path'                     => 0,
        'File::ReadBackwards'            => 0,
        'File::Spec'                     => 0,
        'File::Temp'                     => 0,
        'Filesys::DfPortable'            => 0,
        'HTTP::Parser::XS'               => 0,
        'Inline::C'                      => 0,
        'Inline::Filters'                => 0,
        'IO::Capture::Stderr'            => 0,
        'IO::Uncompress::AnyUncompress'  => 0,
        'List::Util'                     => 0,
        'List::MoreUtils'                => 0,
        'LWP::UserAgent'                 => 0,
        'Module::Find'                   => 0,
        'Moose'                          => 0,
        'MooseX::Aliases'                => 0,
        'MooseX::Declare'                => 0,
        'MooseX::NonMoose'               => 0,
        'MooseX::StrictConstructor'      => 0,
        'MooseX::Types'                  => 0,
        'MooseX::Types::Parameterizable' => '>= 0.07, != 0.08', # 0.08, currently the latest version, just doesn't seem to pass tests for anyone
        'Net::FTP::Robust'               => 0,
        'Net::SSH'                       => 0,
        'Parallel::ForkManager'          => 0,
        'Path::Class'                    => 0,
        'Perl6::Form'                    => 0,
        'POSIX'                          => 0,
        'Proc::Killfam'                  => 0,
        'Storable'                       => 0,
        'Sys::CPU'                       => 0,
        'Sys::Hostname'                  => 0,
        'Sys::Hostname::Long'            => 0,
        'Test::DBIx::Class'              => 0,
        'Test::Most'                     => 0,
        'Test::Strict'                   => 0,
        'Time::Format'                   => 0,
        'Twiggy'                         => 0,
        'TryCatch'                       => 0,
        'XML::LibXML'                    => 0,
        'XML::LibXSLT'                   => 0
    };
}

our %do_not_use = (
    'perl'                           => 1,
    'DBIx::Class::DeploymentHandler' => 1,
    'Inline::C'                      => 1
);

sub create_site_config {
    my $self = shift;
    
    # VRPipe::Config has CPAN dependencies, but we want Module::Build to let the
    # user auto-install those, so we can't just use V:C normally
    my $error;
    foreach my $module (keys %{ required_modules() }, 'VRPipe::Config') {
        next if exists $do_not_use{$module};
        next if $module =~ /^Test/;
        eval "use $module;";
        if ($@) {
            $error = "The actual error when trying to use $module:\n" . $@;
            last;
        }
    }
    if ($error) {
        $error = "

WARNING WARNING WARNING
There is a problem creating the SiteConfig - you're probably missing some
essential dependency modules. Be sure to run ./Build installdeps and then run
this script again, or nothing will work!
WARNING WARNING WARNING

$error
";
        warn $error;
        return;
    }
    
    my $vrp_config             = VRPipe::Config->new();
    my $siteconfig_module_path = $vrp_config->config_module_path;
    
    my $do_config = 'y';
    if (-s $siteconfig_module_path) {
        $do_config = $self->y_n('SiteConfig module already exists; do you wish to go through setup again? [y|n]', 'n');
    }
    
    if ($do_config) {
        print "
When answering the following questions, hit return to accept the default given in [].
To say that you don't want the default or any other value, type the word 'undef' without quotes.
To specify the answer should come from an environment variable, type 'ENV{variable_name}' without the quotes, replacing 'variable_name' as desired\n\n";
        
        while (my $option = $vrp_config->next_option) {
            $option->prompt;
        }
        
        $vrp_config->write_config_module();
        warn "Wrote configuration to $siteconfig_module_path\n";
    }
}

sub ACTION_realclean {
    my $self = shift;
    
    eval 'use Path::Class;';
    eval 'use VRPipe::Config;';
    unless ($@) {
        my $vrp_config             = VRPipe::Config->new();
        my $siteconfig_module_path = $vrp_config->config_module_path;
        
        if (dir()->absolute->contains($siteconfig_module_path->dir)) {
            unlink($siteconfig_module_path);
        }
    }
    
    $self->SUPER::ACTION_realclean(@_);
}

sub ACTION_test {
    my $self = shift;
    
    # we must start vrpipe-server before testing:
    my $local_script = File::Spec->catfile('scripts', 'vrpipe-server');
    my $server;
    if (-x $local_script && -d 'modules' && -d 't') {
        my $thisperl = $Config{perlpath};
        if ($^O ne 'VMS') {
            $thisperl .= $Config{_exe} unless $thisperl =~ m/$Config{_exe}$/i;
        }
        $server = "$thisperl -Imodules -It $local_script --deployment testing";
    }
    unless ($server) {
        die "Can't run tests because this doesn't seem to be the root of the git repository (executable scripts and/or modules directory were not found)\n";
    }
    
    system("$server --farm testing_farm restart");
    warn "If tests are interrupted, you can manually stop the server with this command: $server stop\n";
    
    eval { $self->SUPER::ACTION_test(@_); };
    
    # and then stop it afterwards. This is the real point of this override:
    # the testing server might stay alive on the users system until their next
    # reboot, which is wasteful if they never do any development or testing!
    system("$server stop");
}

# this method, called by create_build_script(), is just annoying and we don't
# need what it does, at least for now
sub create_mymeta {
    return 1;
}

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
        $thing = $dir . "/" . $thing;
        
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

1;
