=head1 NAME

VRPipe::Base::Configuration::Trait::Attribute::ConfigKey - config keys

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

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

use VRPipe::Base;

role VRPipe::Base::Configuration::Trait::Attribute::ConfigKey {
    has question => (
        is        => 'ro',
        isa       => 'Str',
        predicate => 'has_question'
    );
    
    has question_number => (
        is        => 'ro',
        isa       => 'Str',
        required  => 1
    );
    
    has valid => (
        is        => 'ro',
        isa       => ArrayRefOfStrings,
        predicate => 'has_valid'
    );
    
    has env => (
        is        => 'ro',
        isa       => 'Str',
        predicate => 'has_env'
    );
    
    has skip => (
        is        => 'ro',
        isa       => 'Str',
        predicate => 'has_skip'
    );
    
    has secure => (
        is        => 'ro',
        isa       => 'Bool',
        predicate => 'has_secure'
    );
    
    around _process_options (ClassName|Object $class: Str $name, HashRef $options) {
        $options->{isa} = MaybeStrOrEnv;
        $class->$orig($name, $options);
        $options->{lazy} = 1;
        $class->_process_default_or_builder_option($name, $options);
    }
    
    sub _process_default_or_builder_option {
        my $class   = shift;
        my $name    = shift;
        my $options = shift;
        
        my $env_key = $options->{env};
        unless ($env_key) {
            $env_key = 'vrpipe_'.$name;
        }
        
        if (exists $options->{default}) {
            my $def = $options->{default};
            
            if (ref $options->{default}) {
                $options->{default} = sub {
                    my $val = $_[0]->_from_config_or_env($name, $env_key);
                    return $val if defined $val;
                    return $def->($_[0]);
                };
            }
            else {
                $options->{default} = sub {
                    my $val = $_[0]->_from_config_or_env($name, $env_key);
                    return $val if defined $val;
                    return $def;
                };
            }
        }
        elsif ($options->{builder}) {
            my $builder = delete $options->{builder};
            $options->{default} = sub {
                my $val = $_[0]->_from_config_or_env($name, $env_key);
                return $val if defined $val;
                return $_[0]->$builder();
            };
        }
        else {
            $options->{default} = sub {
                return $_[0]->_from_config_or_env($name, $env_key);
            };
        }
    }
}

1;