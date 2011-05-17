package VRPipe::Base::Moose;
use Moose;

extends 'Moose::Object';
with 'VRPipe::Base::Debuggable';
with 'VRPipe::Base::Cleanable';
with 'VRPipe::Base::FileMethods';

no Moose;

1;