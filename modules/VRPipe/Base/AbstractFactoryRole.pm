# stolen directly from MooseX::AbstractFactory::Role, with tweaks to avoid
# memory leak
#
# This file is copyright (c) 2012 by Mike Whitaker.
#
# This is free software; you can redistribute it and/or modify it under the same
# terms as the Perl 5 programming language system itself.

package VRPipe::Base::AbstractFactoryRole;
use strict;
use warnings;
use Moose::Role;
use TryCatch;

has _options        => (is => 'ro', isa => 'ArrayRef[Any]');
has _implementation => (is => 'ro', isa => 'Str');

sub create {
    my ($class, $impl, @impl_args) = @_;
    
    if (defined $impl) {
        my $factory = $class->new({
                _implementation => $impl,
                _options        => [@impl_args]
            }
        );
        
        my $iclass = $factory->_get_implementation_class($factory->_implementation());
        
        # pull in our implementation class
        $factory->_validate_implementation_class($iclass);
        
        my $iconstructor = $iclass->meta->constructor_name;
        
        my $implementation = $iclass->$iconstructor(@{ $factory->_options });
        
        return $implementation;
    }
    else {
        confess('No implementation provided');
    }
}

sub _get_implementation_class {
    my ($self, $impl) = @_;
    
    my $class = blessed $self;
    if ($self->meta->has_class_maker) {
        return $self->meta->implementation_class_maker->($impl);
    }
    else {
        return $class . "::$impl";
    }
}

sub _validate_implementation_class {
    my ($self, $iclass) = @_;
    
    my $err;
    try {
        # can we load the class?
        eval "use $iclass;"; # may die if user really stuffed up _get_implementation_class()
        die $@ if $@;
        
        # does it do the correct roles?
        if ($self->meta->has_implementation_roles) {
            foreach my $role (@{ $self->meta->implementation_roles() }) {
                die "does not implement role $role\n" unless $iclass->does($role);
            }
        }
    }
    catch ($err) {
        confess "Invalid implementation class $iclass: $err";
    }
    
    return;
}

1;
