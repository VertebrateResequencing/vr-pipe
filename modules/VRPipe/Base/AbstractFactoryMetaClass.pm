# stolen directly from MooseX::AbstractFactory::Meta::Class

package VRPipe::Base::AbstractFactoryMetaClass;
use strict;
use warnings;
use Moose;
extends 'Moose::Meta::Class';

has implementation_roles => (
    isa       => 'ArrayRef',
    is        => 'rw',
    predicate => 'has_implementation_roles',
);

has implementation_class_maker => (
    isa       => 'CodeRef',
    is        => 'rw',
    predicate => 'has_class_maker',
);

1;
