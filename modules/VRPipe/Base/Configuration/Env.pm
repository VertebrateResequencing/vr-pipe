use VRPipe::Base;

class VRPipe::Base::Configuration::Env {
    use overload(q[""] => 'stringify', fallback => 1);
    __PACKAGE__->meta->add_package_symbol('&()' => sub { });
    __PACKAGE__->meta->add_package_symbol('&(""' => sub { shift->stringify });
    
    has variable => (
        is        => 'ro',
        isa       => 'Str',
        required  => 1
    );
    
    # we don't have a value attribute so that the value is never stored in the
    # config file
    
    method value {
        return $ENV{$self->variable};
    }
    
    method stringify {
        $self->value ? $self->value : '';
    }
}