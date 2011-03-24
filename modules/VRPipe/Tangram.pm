=head1 DESCRIPTION

Just testing out Tangram

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

package VRPipe::Tangram;

use Heritable::Types;
use Tangram::Core;
use Tangram::Type::Dump::Any;

our $schema =
    Tangram::Schema->new( { classes =>
                            [ HASH => {
                                    fields => {
                                        idbif => undef
                                    }
                                }
                            ]
                           } );

sub db {
    Tangram::Storage->new($schema, @_);
}

1;
