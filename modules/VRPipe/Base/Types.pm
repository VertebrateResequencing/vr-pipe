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
                                  ArrayRefOfStrings FileOrHandle Varchar
                                  IntSQL File Dir MaybeFile MaybeDir StrOrEnv)];

# import built-in types to subtype from
use MooseX::Types::Parameterizable qw(Parameterizable);
use MooseX::Types::Moose qw(Any Num Int Defined Str FileHandle ArrayRef HashRef Maybe);
use Path::Class;

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

class_type('VRPipe::Base::Configuration::Env');
subtype StrOrEnv,
    as 'Str | VRPipe::Base::Configuration::Env',
    #where { Str->check($_) || 'VRPipe::Base::Configuration::Env'->check($_) },
    message { "$_ is neither a String nor a VRPipe::Base::Configuration::Env" };

# file-related (mostly stolen from MooseX::Types::Path::Class)
class_type('Path::Class::Dir');
class_type('Path::Class::File');

subtype Dir, as 'Path::Class::Dir', where { "$_" =~ /^[-\w.#\/\\]+$/ }, message { defined $_ ? "'$_' does not seem like a directory" : "no directory specified" };
subtype File, as 'Path::Class::File', where { "$_" =~ /^[-\w.#\/\\]+$/ }, message { defined $_ ? "'$_' does not seem like a file" : "no file specified" };
subtype MaybeFile, as Maybe[File];
subtype MaybeDir, as Maybe[Dir];
    
for my $type ('Path::Class::Dir', Dir, MaybeDir) {
    coerce $type,
        from Str,      via { Path::Class::Dir->new($_) },
        from ArrayRef, via { Path::Class::Dir->new(@{$_}) };
}

for my $type ('Path::Class::File', File, MaybeFile) {
    coerce $type,
        from Str,      via { Path::Class::File->new($_) },
        from ArrayRef, via { Path::Class::File->new(@{$_}) };
}

subtype FileOrHandle,
    as Defined,
    where { File->check($_) || FileHandle->check($_) },
    message { "Neither a file name nor handle" };

# Varchar and IntSQL for database constraints
subtype Varchar,
    as Parameterizable[Str, Int],
    where {
        my ($string, $int) = @_;
        $int >= length($string) ? 1 : 0;
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
