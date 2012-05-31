=head1 NAME

VRPipe::Base::Configuration - provides configuration mechanics

=head1 SYNOPSIS

    use VRPipe::Base;
    
    class VRPipe::Config {
        use VRPipe::Base::Configuration;
        
        my $question_number = 0;
        
        has dbtype => (
            is      => 'rw',
            question => 'What DRM should be used?',
            default => 'mysql',
            env => 'DBI_DRIVER',
            valid => [qw(mysql postgres)],
            question_number => ++$question_number
        );
        
        has dbname => (
            is      => 'rw',
            question => 'What is the name of your database?',
            question_number => ++$question_number
        );
        
        has db_username => (
            is      => 'rw',
            question => 'What username is used to connect to your database?',
            default => '',
            env     => 'DB_USER',
            question_number => ++$question_number
        );
        
        has db_password => (
            is      => 'rw',
            question => 'What password is used to connect to your database?',
            default => '',
            env     => 'DB_PASSWORD',
            secure => 1, # password will be stored encrypted
            question_number => ++$question_number
        );
    }

    package main;
    my $config = VRPipe::Config->new();
    
    my $dbtype= $config->dbtype;
    $config->dbname('foo');
    
    # loop through all the options (dbtype, dbname, db_username, db_password):
    while (my $option = $config->next_option) {
        # $option isa VRPipe::Base::Configuration::Option
        my $question = $option->question;
        my $question_number = $option->question_number;
        my $option_key = $option->key;
        
        # get/set an option
        my $value = $option->value('new_value');
        
        # prompt user on the command line for the answer to the question
        $option->prompt;
    }
    
    $config->write_config_module();

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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
