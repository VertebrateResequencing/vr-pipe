=pod

http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/ResultSource.pm#add_columns

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
}

1;