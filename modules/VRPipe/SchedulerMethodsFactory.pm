package VRPipe::SchedulerMethodsFactory;
use MooseX::AbstractFactory;

implementation_does qw/VRPipe::SchedulerMethodsRole/;
implementation_class_via sub { 'VRPipe::Schedulers::'.shift };

1;