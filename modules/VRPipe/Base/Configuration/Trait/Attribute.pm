use VRPipe::Base;

role VRPipe::Base::Configuration::Trait::Attribute {
    around interpolate_class ($class: HashRef $options) {
        $options->{traits} ||= [];
        push @{$options->{traits}}, 'VRPipe::Base::Configuration::Trait::Attribute::ConfigKey';
        return $class->$orig($options);
    }
}

1;