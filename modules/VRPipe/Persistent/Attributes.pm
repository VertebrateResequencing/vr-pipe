
=head1 NAME

VRPipe::Persistent::Attributes - a role for Persistent, specifying attributes

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

See:
L<http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/ResultSource.pm#add_columns>.

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

role VRPipe::Persistent::Attributes {
    has is_auto_increment => (
        is        => 'rw',
        isa       => 'Bool',
        default   => 0,
        predicate => 'is_auto_increment_was_set'
    );
    has is_primary_key => (
        is        => 'rw',
        isa       => 'Bool',
        default   => 0,
        predicate => 'is_primary_key_was_set'
    );
    has is_key => (
        is        => 'rw',
        isa       => 'Bool',
        default   => 0,
        predicate => 'is_key_was_set'
    );
    has allow_key_to_default => (
        is        => 'rw',
        isa       => 'Bool',
        default   => 0,
        predicate => 'allow_key_to_default_was_set'
    );
    has _key_default => (
        is        => 'rw',
        isa       => 'Defined',
        predicate => '_key_default_was_set'
    );
    has is_nullable => (
        is        => 'rw',
        isa       => 'Bool',
        default   => 0,
        predicate => 'is_nullable_was_set'
    );
    has default_value => (
        is        => 'rw',
        isa       => 'Defined',
        predicate => 'default_value_was_set'
    );
    has extra => (
        is        => 'rw',
        isa       => 'HashRef',
        predicate => 'extra_was_set'
    );
    
    # relationships
    has [qw(belongs_to has_one might_have)] => (is => 'rw', isa => RelationshipArg);
    
    around _process_options (ClassName|Object $class: Str $name, HashRef $options) {
        $options->{clearer} = '_clear_' . $name;
        my $isa = $options->{isa};
        if ($isa eq 'HashRef' || $isa eq 'ArrayRef') {
            $options->{isa} = "$isa\[ArrayRef[Str]|HashRef[Str]|Str]|Str"; # persistent->get() will freeze the ref before passing to find(), so we must validate as Str as well
        }
        $class->$orig($name, $options);
        $class->_process_default_or_builder_option($name, $options);
    }
    
    sub _process_default_or_builder_option {
        my $class   = shift;
        my $name    = shift;
        my $options = shift;
        
        # before returning default value, set that value in dbic accessor by
        # calling the method accessor which Persistent has a custom around
        # modifier for
        if (exists $options->{default} || exists $options->{builder}) {
            $options->{lazy} = 0; # defaults must not be lazy, or their default values will not get set in the db
            
            if (exists $options->{default}) {
                my $def = $options->{default};
                
                if (ref $options->{default}) {
                    $options->{default} = sub {
                        my $self   = shift;
                        my $return = $def->($self);
                        if (ref($self)) {
                            $self->$name($return);
                        }
                        return $return;
                    };
                }
                # else, the database will have the correct default, so dbic will
                # get the correct value from the db and we don't need to set it
                # via our accessor
            }
            elsif (exists $options->{builder}) {
                my $builder = delete $options->{builder};
                $options->{default} = sub {
                    my $self   = shift;
                    my $return = $self->$builder();
                    if (ref($self)) {
                        $self->$name($return);
                    }
                    return $return;
                };
            }
            
            if ($options->{allow_key_to_default}) {
                $options->{_key_default} = $options->{default};
            }
        }
    }
}

1;
