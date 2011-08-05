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
                     isa => Text,
                     coerce => 1,
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'options' => (is => 'rw',
                      isa => 'HashRef',
                      traits => ['VRPipe::Persistent::Attributes'],
                      default => sub { {} },
                      allow_key_to_default => 1,
                      is_key => 1);
    
    has 'description' => (is => 'rw',
                          isa => Text,
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1);
    
    has '_changed_marker' => (is => 'rw',
                             isa => Varchar[255],
                             traits => ['VRPipe::Persistent::Attributes'],
                             is_nullable => 1);
    
    has '_source_instance' => (is => 'rw',
                               isa => 'Defined',
                               lazy => 1,
                               builder => '_build_source');
    
    has '_elements' => (is => 'rw',
                        isa => 'ArrayRef',
                        lazy => 1,
                        builder => '_build_elements');
    
    has '_next_element' => (is => 'rw',
                            isa => 'Int',
                            default => 0);
    
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
        
        my $next = $self->_next_element || 0;
        
        # if we've never parsed the source, or it has changed, we'll parse it
        # now and get it to generate all our dataelements
        if ($next == 0 && $source->_has_changed) {
            my $latest_elements = $source->_get_elements;
            
            # if we already have elements from a previous parsing of the source,
            # check if any are no longer desired
            my $current_elements = $self->_elements;
            if (@$current_elements) {
                my %good_ids = map { $_->id => 1 } @$latest_elements;
                foreach my $element (@$current_elements) {
                    unless (exists $good_ids{$element->id}) {
                        $element->withdrawn(1);
                        $element->update;
                    }
                }
            }
            $self->_elements($self->_build_elements);
            $self->_changed_marker($source->_changed_marker);
            $self->update;
        }
        
        my $elements = $self->_elements;
        if ($next > $#$elements) {
            return;
        }
        $self->_next_element($next + 1);
        return $elements->[$next];
    }
    
    method _build_source {
        my $changed_marker = $self->_changed_marker;
        return VRPipe::DataSourceFactory->create($self->type, {method => $self->method,
                                                               source => $self->source,
                                                               options => $self->options,
                                                               $changed_marker ? ('_changed_marker' => $changed_marker) : (),
                                                               '_datasource_id' => $self->id});
    }
    
    method _build_elements {
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('DataElement')->search({ datasource => $self->id, withdrawn => 0 });
        
        # we can't store this $rs and call $rs->next in next_element(), because
        # something else might call disconnect, breaking our search.
        # instead we get all elements up-front and store an array ref of
        # them, returning the next one each time. This is actually very fast,
        # even with ~10000 data elements
        
        my @elements;
        while (my $element = $rs->next) {
            push(@elements, $element);
        }
        return \@elements;
    }
}

1;