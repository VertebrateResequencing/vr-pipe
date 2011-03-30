package VRPipe::Base::Moose;
use Moose;

extends 'Moose::Object';
with 'VRPipe::Base::Debuggable';
with 'VRPipe::Base::Cleanable';

sub foo { return 'food'; }

no Moose;

1;