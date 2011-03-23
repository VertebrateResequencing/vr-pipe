=head1 NAME

VRPipe::Base::Types - MooseX type definitions

=head1 SYNOPSIS

package Foo;
use Moose;
use VRPipe::Base::Types qw(PositiveInt VerbosityValue);

# use the exported constants as type names
has 'bar',
    isa    => PositiveInt,
    is     => 'rw';
has 'baz',
    isa    => VerbosityValue,
    is     => 'rw';

=head1 DESCRIPTION

Type definitions for use by VRPipe classes.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

package VRPipe::Base::Types;
use strict;
use warnings;

# predeclare our own types
use MooseX::Types -declare => [qw(PositiveInt VerbosityValue ArrayRefOfInts
                                  ArrayRefOfStrings FileNameOrHandle Varchar
                                  IntSQL)];

# import built-in types to subtype from
use MooseX::Types::Parameterizable qw(Parameterizable);
use MooseX::Types::Moose qw(Num Int Defined Str FileHandle ArrayRef);

# custom type definitions
subtype PositiveInt, 
    as Int, 
    where { $_ > 0 },
    message { "Int is not larger than 0" };
coerce PositiveInt,
    from Int,
    via { 1 };

subtype VerbosityValue,
    as Num,
    where { $_ == -1 || $_ == 0 || $_ == 0.5 || $_ == 1 || $_ == 2 },
    message { "Verbosity is not amongst valid values -1|0|0.5|1|2" };

#subtype File,
#    as Object,
#    where { $_->does('VRPipe::Base::ReadWritable') };
#coerce File,
#    from FileHandle,
#    via { VRPipe::Base::FHWrapper->new(handle => $_) };

subtype FileNameOrHandle,
    as Defined,
    where { Str->check($_) || FileHandle->check($_) },
    message { "Neither a file name nor handle" };

# Varchar and IntSQL for database constraints
subtype Varchar,
    as Parameterizable[Str, Int],
    where {
        my ($string, $int) = @_;
        $string eq 'Footon' ? 1 : 1; #$int >= length($string) ? 1 : 0;
    },
    message { "'$_' is too long" };
subtype IntSQL,
    as Parameterizable[Int, Int],
    where {
        my ($number, $int) = @_;
        (defined $number && $int >= length("$number")) ? 1 : 0;
    },
    message { defined $_ ? "'$_' is too long" : "number is undefined" };


# allow users to supply either a single X, or an array ref of them, eg:
# has 'my_attribute' => ( is => 'rw', isa => 'ArrayRefOfInts', coerce => 1 );
# $obj->new(my_attribute => 42); $obj->new(my_attribute => [42, 31]);
subtype ArrayRefOfInts,
    as ArrayRef[Int];
coerce ArrayRefOfInts,
    from Int,
    via { [ $_ ] };

subtype ArrayRefOfStrings,
    as ArrayRef[Str];
coerce ArrayRefOfStrings,
    from Str,
    via { [ $_ ] };

1;
