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
    
    has '_changed_marker' => (is => 'rw',
                              isa => 'Str');
    
    has '_datasource_id' => (is => 'ro',
                             isa => Persistent);
    
    requires '_open_source';
    requires '_has_changed';
    requires '_update_changed_marker';
    
    method _get_elements {
        my $handle = $self->_handle || return;
        my $method = $self->method;
        $self->can($method) || $self->throw("Invalid method '$method' for ".ref($self));
        my $elements = $self->$method(%{$self->options}, handle => $handle);
        $self->_update_changed_marker;
        return $elements;
    }
}

1;