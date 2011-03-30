package VRPipe::Base::Moose;
use Moose;

extends 'Moose::Object';
with 'VRPipe::Base::Debuggable';
with 'VRPipe::Base::Cleanable';

no Moose;

1;