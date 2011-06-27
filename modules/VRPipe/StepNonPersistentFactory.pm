package VRPipe::StepNonPersistentFactory;
use MooseX::AbstractFactory;

implementation_does qw/VRPipe::StepRole/;
implementation_class_via sub { 'VRPipe::Steps::'.shift };

1;