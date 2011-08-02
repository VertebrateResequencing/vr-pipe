package VRPipe::Base::UseMoose;

use Moose ();
use Moose::Exporter;
use Moose::Util::MetaRole;

Moose::Exporter->setup_import_methods(also => 'Moose');

sub init_meta {
    shift;
    my %args = @_;
    
    Moose->init_meta(%args, base_class => 'VRPipe::Base::Moose');
    
    # for some reason this is required for something like 'after BUILDALL' to
    # work in VRPipe::Base::Debuggable
    #Moose::Util::MetaRole::apply_metaroles(
    #    for             => $args{for_class},
    #    class_metaroles => {
    #        class => => ['MooseX::NonMoose::Meta::Role::Class'],
    #        constructor => ['MooseX::NonMoose::Meta::Role::Constructor'],
    #    },
    #);
    
    # roles can be set like this, or by using 'with' in VRPipe::Base::Moose;
    # not sure which is best?
    #Moose::Util::MetaRole::apply_base_class_roles(
    #    for   => $args{for_class},
    #    roles => ['VRPipe::Base::Debuggable'],
    #);
    
    return $args{for_class}->meta();
}

no Moose;

1;