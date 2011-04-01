=head1 NAME

VRPipe::Base::Configuration - provides configuration mechanics

=head1 SYNOPSIS

use VRPipe::Base;

class VRPipe::Config {
    use VRPipe::Base::Configuration;
    
    has database_name => (
        is      => 'ro',
        isa     => 'Str',
        default => 'MyApp',
        section => 'database',
        key     => 'name',
        documentation =>
            'The name of the database.',
    );
    
    has database_username => (
        is      => 'ro',
        isa     => Str,
        default => q{},
        section => 'database',
        key     => 'username',
        documentation =>
            'The username to use when connecting to the database. By default, this is empty.',
    );
}

package main;
my $config = VRPipe::Config->new();

#...

$config->write_config_module(values => %values);

=head1 DESCRIPTION

Based on MooseX::Configuration
http://search.cpan.org/~drolsky/MooseX-Configuration-0.02/lib/MooseX/Configuration.pm
Modified to read from and write to a Perl module instead of an ini file. Other
interface changes were made; see synopsis.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

package VRPipe::Base::Configuration;

use strict;
use warnings;

use Moose::Exporter;

Moose::Exporter->setup_import_methods(
    class_metaroles => {
        attribute => ['VRPipe::Base::Configuration::Trait::Attribute'],
    },
    base_class_roles => ['VRPipe::Base::Configuration::Trait::Object'],
);

1;
