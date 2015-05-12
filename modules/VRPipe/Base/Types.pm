
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

package VRPipe::Base::Types;
use strict;
use warnings;

# predeclare our own types
use MooseX::Types -declare => [
    qw(PositiveInt VerbosityValue ArrayRefOfInts
      ArrayRefOfStrings FileOrHandle Varchar
      IntSQL File Dir MaybeFile MaybeDir StrOrEnv
      MaybeStrOrEnv Datetime VRPFileOrHandle
      Persistent PersistentObject RelationshipArg
      PersistentArrayRef PersistentHashRef
      FileType FileProtocol AbsoluteFile PersistentFileHashRef
      OpenMode AnyFileHandle
      ParserType MapperType PreviousStepOutput
      Text)
];

# import built-in types to subtype from
use MooseX::Types::Parameterizable qw(Parameterizable);
use MooseX::Types::Moose qw(Any Num Int Defined Str FileHandle ArrayRef HashRef Maybe Object ClassName);
use Path::Class;
use File::HomeDir;
use DateTime ();
use File::Temp;

# custom type definitions
subtype PositiveInt, as Int, where { $_ > 0 }, message { "Int is not larger than 0" };
coerce PositiveInt, from Int, via { 1 };

subtype VerbosityValue, as Num, where { $_ == -1 || $_ == 0 || $_ == 0.5 || $_ == 1 || $_ == 2 }, message { "Verbosity is not amongst valid values -1|0|0.5|1|2" };

class_type('VRPipe::Base::Configuration::Env');
subtype StrOrEnv, as 'Str | VRPipe::Base::Configuration::Env', message { "$_ is neither a String nor a VRPipe::Base::Configuration::Env" };
subtype MaybeStrOrEnv, as Maybe [StrOrEnv];

subtype RelationshipArg, as Defined, where { ClassName->check($_) || ArrayRef->check($_) }, message { "$_ is neither a class name nor an array ref" };

class_type('VRPipe::Persistent::Schema');
class_type('VRPipe::Persistent');
class_type('VRPipe::Submission');
class_type('VRPipe::Requirements');
class_type('VRPipe::StepState');
class_type('VRPipe::File');
class_type('VRPipe::Pipeline');
class_type('VRPipe::PipelineSetup');
class_type('VRPipe::DataElement');
class_type('VRPipe::Step');
class_type('VRPipe::StepCmdSummary');
class_type('VRPipe::StepAdaptorDefiner');
class_type('VRPipe::StepBehaviourDefiner');
class_type('VRPipe::StepMember');
class_type('VRPipe::Interface::BackEnd');
class_type('VRPipe::DataElementState');

# file-related (mostly stolen from MooseX::Types::Path::Class)
class_type('Path::Class::Dir');
class_type('Path::Class::File');

subtype Dir,  as 'Path::Class::Dir',  where { "$_" =~ /^[-+\w.,#\/\\~:]+$/ }, message { defined $_ ? "'$_' does not seem like a directory" : "no directory specified" };
subtype File, as 'Path::Class::File', where { "$_" =~ /^[-+\w.,#\/\\~:]+$/ }, message { defined $_ ? "'$_' does not seem like a file"      : "no file specified" };
subtype MaybeFile, as Maybe [File];
subtype MaybeDir,  as Maybe [Dir];
subtype AbsoluteFile, as File, where { $_->is_absolute }, message {
    my $file = $_;
    if   (ref($file) && $file->isa('Path::Class::File')) { return "path '$_' must be absolute"; }
    else                                                 { return "Not a Path::Class::File"; }
};

sub cleanup_path_str {
    my $str = shift;
    return unless $str;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    if ($str =~ /^~/) {
        my $home = File::HomeDir->my_home;
        $str =~ s/^~/$home/;
    }
    return $str;
}

for my $type ('Path::Class::Dir', Dir, MaybeDir) {
    coerce $type, from Str, via {
        Path::Class::Dir->new(cleanup_path_str($_));
    }, from ArrayRef, via {
        Path::Class::Dir->new(map { cleanup_path_str($_) } @$_);
    };
}

for my $type ('Path::Class::File', File, MaybeFile, AbsoluteFile) {
    coerce $type, from Str, via {
        Path::Class::File->new(cleanup_path_str($_));
    }, from ArrayRef, via { #*** this 'works', but also manages to store a non-functional ARRAY... string if used to create a VRPipe::File
        Path::Class::File->new(map { cleanup_path_str($_) } @$_);
    };
}

class_type('IO::File');
subtype FileOrHandle, as Defined, where { File->check($_) || FileHandle->check($_) || (Object->check($_) && $_->isa('IO::File')) }, message { "Neither a file name nor handle" };
coerce FileOrHandle, from Str, via {
    if (/^~/) { my $home = File::HomeDir->my_home; $_ =~ s/^~/$home/; }
    Path::Class::File->new($_);
};

subtype VRPFileOrHandle, as Defined, where { File->check($_) || FileHandle->check($_) || (Object->check($_) && ($_->isa('IO::File') || $_->isa('VRPipe::File'))) }, message { "Not a file name, handle or VRPipe::File" };
coerce VRPFileOrHandle, from Str, via {
    if (/^~/) { my $home = File::HomeDir->my_home; $_ =~ s/^~/$home/; }
    Path::Class::File->new($_);
};

subtype AnyFileHandle, as Defined, where { FileHandle->check($_) || (Object->check($_) && $_->isa('IO::File')) }, message { "Doesn't look like a file handle" };

class_type('File::Temp::Dir');
class_type('File::Temp::File');

subtype FileType, as Str, where {
    my $type = $_;
    length($type) <= 4 || return 0;
    
    # we used to more do checking, but now anything can be a valid filetype
    # because arbitrary types are allowed
    return 1;
}, message { "Not a valid VRPipe::FileType type" };
coerce FileType, from Str, via { lc($_) };

subtype FileProtocol, as Str, where {
    my $type = $_;
    eval "require VRPipe::FileProtocol::$type;";
    return $@ ? 0 : 1;
}, message { "Not a valid VRPipe::FileProtocol type" };
coerce FileProtocol, from Str, via { lc($_) };

subtype ParserType, as Str, where {
    my $type = $_;
    eval "require VRPipe::Parser::$type;";
    if ($@) { return 0; }
    return 1;
}, message { "Not a valid VRPipe::Parser type" };
coerce ParserType, from Str, via { lc($_) };

subtype MapperType, as Str, where {
    my $type = $_;
    eval "require VRPipe::Wrapper::Mapper::$type;";
    if ($@) { return 0; }
    return 1;
}, message { "Not a valid VRPipe::Wrapper::Mapper type" };
coerce MapperType, from Str, via { lc($_) };

subtype OpenMode, as Str, where { /^(?:>>?|<)$/ }, message { "open modes are restricted to just >, >> and <" };

# datetime (stolen from MooseX::Types::DateTime)
class_type 'DateTime';
subtype Datetime, as 'DateTime';
for my $type ("DateTime", Datetime) {
    coerce $type => from Num, via { 'DateTime'->from_epoch(epoch => $_) }, from HashRef, via { 'DateTime'->new(%$_) };
}

# database constraints
subtype Varchar, as Parameterizable [Str, Int], where {
    my ($string, $int) = @_;
    $int >= length($string) ? 1 : 0;
}, message { "'$_' is too long" };

subtype Text, as Str;
coerce Text, from 'Path::Class::File', via { $_->absolute->stringify };

subtype IntSQL, as Parameterizable [Int, Int], where {
    my ($number, $int) = @_;
    (defined $number && $int >= length("$number")) ? 1 : 0;
}, message { defined $_ ? "'$_' is too long" : "number is undefined" };

subtype PersistentObject, as 'VRPipe::Persistent';
#    as Object,
#    where { $_->isa('VRPipe::Persistent') },
#    message { "Not a Persistent object" };

subtype Persistent, as PositiveInt, # can't coerce IntSQL[16]
  where { length(shift) <= 16 };
coerce Persistent, from PersistentObject, via { $_->{_column_data}->{id} }; # this is called so many times, its worth directly accessing the hash structure instead of calling id(), for the speedup

subtype PersistentArrayRef, as ArrayRef [PersistentObject];

subtype PersistentHashRef, as HashRef [PersistentObject];

subtype PersistentFileHashRef, as HashRef [ArrayRef ['VRPipe::File']];

subtype PreviousStepOutput, as HashRef [HashRef [ArrayRef ['VRPipe::File']]];

# allow users to supply either a single X, or an array ref of them, eg:
# has 'my_attribute' => ( is => 'rw', isa => 'ArrayRefOfInts', coerce => 1 );
# $obj->new(my_attribute => 42); $obj->new(my_attribute => [42, 31]);
subtype ArrayRefOfInts, as ArrayRef [Int];
coerce ArrayRefOfInts, from Int, via { [$_] };

subtype ArrayRefOfStrings, as ArrayRef [Str];
coerce ArrayRefOfStrings, from Str, via { [$_] };

1;
