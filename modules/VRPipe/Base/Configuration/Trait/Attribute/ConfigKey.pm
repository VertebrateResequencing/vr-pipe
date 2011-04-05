use VRPipe::Base;

role VRPipe::Base::Configuration::Trait::Attribute::ConfigKey {
    has question => (
        is        => 'ro',
        isa       => 'Str',
        predicate => 'has_question'
    );
    
    has question_number => (
        is        => 'ro',
        isa       => 'Str',
        required  => 1
    );
    
    has valid => (
        is        => 'ro',
        isa       => ArrayRefOfStrings,
        predicate => 'has_valid'
    );
    
    has env => (
        is        => 'ro',
        isa       => 'Str',
        predicate => 'has_env'
    );
    
    around _process_options (ClassName|Object $class: Str $name, HashRef $options) {
        $options->{isa} = MaybeStrOrEnv;
        $class->$orig($name, $options);
        $options->{lazy} = 1;
        $class->_process_default_or_builder_option($name, $options);
    }
    
    sub _process_default_or_builder_option {
        my $class   = shift;
        my $name    = shift;
        my $options = shift;
        
        my $env_key = $options->{env};
        unless ($env_key) {
            $env_key = 'vrpipe_'.$name;
        }
        
        if (exists $options->{default}) {
            my $def = $options->{default};
            
            if (ref $options->{default}) {
                $options->{default} = sub {
                    my $val = $_[0]->_from_config_or_env($name, $env_key);
                    return $val if defined $val;
                    return $def->( $_[0] );
                };
            }
            else {
                $options->{default} = sub {
                    my $val = $_[0]->_from_config_or_env($name, $env_key);
                    return $val if defined $val;
                    return $def;
                };
            }
        }
        elsif ($options->{builder}) {
            my $builder = delete $options->{builder};
            
            $options->{default} = sub {
                my $val = $_[0]->_from_config_or_env($name, $env_key);
                return $val if defined $val;
                return $_[0]->$builder();
            };
        }
        else {
            $options->{default} = sub {
                return $_[0]->_from_config_or_env($name, $env_key);
            };
        }
    }
}

1;