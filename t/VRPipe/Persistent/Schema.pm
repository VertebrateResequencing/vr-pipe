use VRPipe::Base;

class t::VRPipe::Persistent::Schema extends DBIx::Class::Schema {
    __PACKAGE__->load_classes({'t::VRPipe' => [qw/Artist CD Track/]});
    #__PACKAGE__->exception_action(sub { die "moo(@_)\n" });
    __PACKAGE__->stacktrace(1);
}

1;
