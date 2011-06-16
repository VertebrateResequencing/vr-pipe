package VRPipe::DataSourceFactory;
use MooseX::AbstractFactory;

implementation_does qw/VRPipe::DataSourceRole/;
implementation_class_via sub { 'VRPipe::DataSource::'.shift };

1;