use VRPipe::Base;

class VRPipe::Base::Configuration::Option {
    use VRPipe::Base::Configuration::Env;
    
    has _attr => (
        is        => 'ro',
        does      => 'VRPipe::Base::Configuration::Trait::Attribute::ConfigKey',
        required  => 1
    );
    
    has _obj => (
        is        => 'ro',
        does      => 'VRPipe::Base::Configuration::Trait::Object',
        weak_ref  => 1,
        required  => 1
    );
    
    method key {
        $self->_attr->name;
    }
    
    method value (StrOrEnv $set?) {
        my $method_name = $self->key;
        
        if ($set) {
            if ($set eq 'undef') {
                return $self->_obj->$method_name(undef);
            }
            
            my $valid = $self->valid;
            if ($valid) {
                my %allowed = map { $_ => 1 } @{$valid};
                if (! exists $allowed{$set}) {
                    $self->throw("'$set' was not a valid value for option '$method_name'");
                }
            }
            return $self->_obj->$method_name($set);
        }
        
        $self->_obj->$method_name();
    }
    
    method env (Str $var) {
        my $env = VRPipe::Base::Configuration::Env->new(variable => $var);
        $self->value($env);
        return $env;
    }
    
    method question {
        my $attr = $self->_attr;
        if ($attr->has_question) {
            return $attr->question;
        }
        else {
            my $key = $self->key;
            return "$key?";
        }
    }
    
    method question_number {
        $self->_attr->question_number;
    }
    
    method valid {
        my $attr = $self->_attr;
        return $attr->valid if $attr->has_valid;
        return;
    }
}

1;