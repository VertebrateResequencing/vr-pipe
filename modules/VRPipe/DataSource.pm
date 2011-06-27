use VRPipe::Base;

class VRPipe::DataSource extends VRPipe::Persistent {
    use VRPipe::DataSourceFactory;
    
    has 'type' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'method' => (is => 'rw',
                     isa => Varchar[64],
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'source' => (is => 'rw',
                     isa => Varchar[64],
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'options' => (is => 'rw',
                      isa => 'HashRef',
                      traits => ['VRPipe::Persistent::Attributes'],
                      default => sub { {} },
                      allow_key_to_default => 1,
                      is_key => 1);
    
    has 'description' => (is => 'rw',
                         isa => Varchar[256],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    has '_source_instance' => (is => 'rw',
                               isa => 'Defined',
                               lazy => 1,
                               builder => '_build_source');
    
    __PACKAGE__->make_persistent(has_many => [elements => 'VRPipe::DataElement']);
    
    around elements {
        # we don't return $self->$orig because that returns all associated
        # elements; we need to first create all elements, and then only return
        # elements that are still "current" (haven't been deleted in the source)
        my @elements;
        while (my $element = $self->next_element) {
            push(@elements, $element);
        }
        return @elements;
    }
    
    method next_element {
        my $source = $self->_source_instance || return;
        my $result = $source->next_result || return;
        return VRPipe::DataElement->get(datasource => $self, result => $result);
    }
    
    method _build_source {
        return VRPipe::DataSourceFactory->create($self->type, {method => $self->method,
                                                               source => $self->source,
                                                               options => $self->options});
    }
}

1;