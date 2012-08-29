# stolen directly from MooseX::AbstractFactory::Meta::Class
#
# This file is copyright (c) 2012 by Mike Whitaker.
#
# This is free software; you can redistribute it and/or modify it under the same
# terms as the Perl 5 programming language system itself.

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
