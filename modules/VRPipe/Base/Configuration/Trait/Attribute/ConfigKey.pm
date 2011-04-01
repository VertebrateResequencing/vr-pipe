use VRPipe::Base;

role VRPipe::Base::Configuration::Trait::Attribute::ConfigKey {
    around _process_options (ClassName|Object $class: Str $name, HashRef $options) {
        $class->$orig($name, $options);
        $class->_process_default_or_builder_option($name, $options);
    }
    
    sub _process_default_or_builder_option {
        my $class   = shift;
        my $name    = shift;
        my $options = shift;
        
        $options->{lazy} = 1;
        
        if (exists $options->{default}) {
            my $def = $options->{default};
            
            if (ref $options->{default}) {
                $options->{default} = sub {
                    my $val = $_[0]->_from_config($name);
                    return $val if defined $val;
                    return $def->( $_[0] );
                };
            }
            else {
                $options->{default} = sub {
                    my $val = $_[0]->_from_config($name);
                    return $val if defined $val;
                    return $def;
                };
            }
        }
        elsif ($options->{builder}) {
            my $builder = delete $options->{builder};
            
            $options->{default} = sub {
                my $val = $_[0]->_from_config($name);
                return $val if defined $val;
                return $_[0]->$builder();
            };
        }
        else {
            $options->{default} = sub {
                return $_[0]->_from_config($name);
            };
        }
    }
}

1;