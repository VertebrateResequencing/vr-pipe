use VRPipe::Base;

role VRPipe::DataSourceRole {
    has 'method' => (is => 'ro',
                     isa => 'Str');
    
    has 'source' => (is => 'ro',
                     isa => 'Defined');
    
    has 'options' => (is => 'ro',
                      isa => 'HashRef');
    
    has '_handle' => (is => 'rw',
                      isa => 'Defined',
                      lazy => 1,
                      builder => '_open_source');
    
    has '_changed_marker' => (is => 'rw',
                              isa => 'Str');
    
    has '_datasource_id' => (is => 'ro',
                             isa => Persistent);
    
    requires '_open_source';
    requires '_has_changed';
    requires '_update_changed_marker';
    requires 'description';
    requires 'source_description';
    requires 'method_description';
    
    method _get_elements {
        my $handle = $self->_handle || return;
        my $method = $self->method;
        $self->can($method) || $self->throw("Invalid method '$method' for ".ref($self));
        my $elements = $self->$method(%{$self->options}, handle => $handle);
        $self->_update_changed_marker;
        return $elements;
    }
    
    method get_methods {
        my $metarole = Moose::Meta::Role->initialize(__PACKAGE__);
        my @meta_methods = $metarole->get_method_list;
        push(@meta_methods, $metarole->get_attribute_list);
        push(@meta_methods, $metarole->get_required_method_list);
        my %role_methods = map { $_ => 1 } @meta_methods;
        
        my $class = ref($self);
        my $classmeta = $self->meta;
        my %self_methods = map { $_->name => 1 } $classmeta->get_all_attributes; #*** really we want to skip methods on all of the roles a class uses, but I don't know how to get that list of all roles...
        $self_methods{new} = 1;
        $self_methods{DESTROY} = 1; #*** bleugh, is there a way to detect these kinds of things automatically as well?
        return map { $_->name } grep { $_->original_package_name eq $class && ! exists $role_methods{$_->name} && ! exists $self_methods{$_->name} && index($_->name, '_') != 0 } $classmeta->get_all_methods;
    }
    
    method method_options (Str $method) {
        $self->throw("'$method' is not a method of ".ref($self)) unless $self->can($method);
        my $method_meta = $self->meta->find_method_by_name($method);
        $self->throw("method $method in class ".ref($self)." is not a signature method") unless $method_meta->isa('MooseX::Method::Signatures::Meta::Method');
        my $sig = $method_meta->parsed_signature;
        
        my @return;
        foreach my $kind (qw(positional named)) {
            my $check_method = 'has_'.$kind.'_params';
            next unless $sig->$check_method;
            # *** actually, we don't support positional args in method options?
            next if $kind eq 'positional';
            my $sig_method = $kind.'_params';
            foreach my $param ($sig->$sig_method) {
                my $var = $param->variable_name;
                next if $var eq '$handle';
                $var =~ s/^\$//; #*** we only handle $ sigils atm...
                my $tc_name = $param->meta_type_constraint->name;
                $tc_name =~ s/VRPipe::Base::Types:://g;
                push(@return, [$kind, $var, $param->required, $param->default_value, $tc_name]);
            }
        }
        
        return @return;
    }
}

1;