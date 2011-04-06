use VRPipe::Base;

class t::VRPipe::Persistent::Schema extends VRPipe::Persistent::SchemaBase {
    __PACKAGE__->load_classes({'t::VRPipe' => [qw/Artist CD Track/]});
}

1;
