=head1 DESCRIPTION

Schema for DBIx::Class

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use VRPipe::Base;

class VRPipe::Persistent::Schema extends DBIx::Class::Schema {
    #__PACKAGE__->load_namespaces();
    #__PACKAGE__->load_namespaces(result_namespace => ['+VRPipe::DirA', '+VRPipe::DirB']);
    #__PACKAGE__->load_classes({VRPipe => [qw/Artist CD Track/]});
    
    # __PACKAGE__->exception_action(sub { My::ExceptionClass->throw(@_) });
    __PACKAGE__->stacktrace(1);
}

1;
