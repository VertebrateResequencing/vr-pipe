package VRPipe::PipelineNonPersistentFactory;
use MooseX::AbstractFactory;

implementation_does qw/VRPipe::PipelineRole/;
implementation_class_via sub { 'VRPipe::Pipelines::'.shift };

1;