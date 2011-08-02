use VRPipe::Base;

role VRPipe::DataSourceRole {
    has 'method' => (is => 'ro',
                     isa => 'Str',
                     required => 1); #*** each method needs a way to self-document itself so front-ends know what it does...
    
    has 'source' => (is => 'ro',
                     isa => 'Defined',
                     required => 1);
    
    has 'options' => (is => 'ro',
                      isa => 'HashRef',
                      required => 1);
    
    has '_handle' => (is => 'rw',
                      isa => 'Defined',
                      lazy => 1,
                      builder => '_open_source');
    
    requires '_open_source';
    
    method next_result {
        my $handle = $self->_handle || return;
        my $method = $self->method;
        $self->can($method) || $self->throw("Invalid method '$method' for ".ref($self));
        return $self->$method(%{$self->options}, handle => $handle);
    }
}

1;