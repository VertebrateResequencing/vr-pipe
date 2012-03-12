package VRPipe::Persistent::ConverterFactory;
use MooseX::AbstractFactory;

implementation_does qw/VRPipe::Persistent::ConverterRole/;
implementation_class_via sub { 'VRPipe::Persistent::Converter::'.shift };

1;
