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
    
    method _build_source {
        my $changed_marker = $self->_changed_marker;
        return VRPipe::DataSourceFactory->create($self->type, {method => $self->method,
                                                               source => $self->source,
                                                               options => $self->options,
                                                               $changed_marker ? ('_changed_marker' => $changed_marker) : (),
                                                               '_datasource_id' => $self->id});
    }
    
    __PACKAGE__->make_persistent(has_many => [elements => 'VRPipe::DataElement']);
    
    around elements {
        # we don't return $self->$orig because that returns all associated
        # elements; we need to first create all elements, and then only return
        # elements that are still "current" (haven't been deleted in the source)
        $self->_prepare_elements_and_states || return;
        
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('DataElement')->search({ datasource => $self->id, withdrawn => 0 });
        
        my @elements;
        while (my $element = $rs->next) {
            push(@elements, $element);
        }
        
        return \@elements;
    }
    
    method incomplete_element_states (VRPipe::PipelineSetup $setup, Int $limit?) {
        $self->_prepare_elements_and_states || return;
        
        my $pipeline = $setup->pipeline;
        my $num_steps = $pipeline->steps;
        
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('DataElementState')->search({ pipelinesetup => $setup->id, completed_steps => {'<', $num_steps}, 'dataelement.withdrawn' => 0 },
                                                                { join => 'dataelement', $limit ? (rows => $limit) : () });
        
        my @incomplete;
        while (my $state = $rs->next) {
            push(@incomplete, $state);
        }
        
        return \@incomplete;
    }
    
    method _prepare_elements_and_states {
        my $source = $self->_source_instance || return;
        
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('PipelineSetup')->search({ datasource => $self->id });
        my @setups;
        while (my $ps = $rs->next) {
            push(@setups, $ps);
        }
        
        if ($source->_has_changed) {
            my $elements = $source->_get_elements;
            
            my %current_elements;
            foreach my $element (@$elements) {
                $current_elements{$element->id} = 1;
                foreach my $setup (@setups) {
                    VRPipe::DataElementState->get(pipelinesetup => $setup, dataelement => $element);
                }
            }
            
            # withdrawn any elements that are no longer in the datasource
            my $schema = $self->result_source->schema;
            my $rs = $schema->resultset('DataElement')->search({ datasource => $self->id, withdrawn => 0 });
            while (my $element = $rs->next) {
                unless (exists $current_elements{$element->id}) {
                    $element->withdrawn(1);
                    $element->update;
                }
            }
            
            $self->_changed_marker($source->_changed_marker);
            $self->update;
        }
        else {
            # check that a new pipelinesetup hasn't been created since the source
            # last changed
            my $expected_count = $schema->resultset('DataElement')->count({ datasource => $self->id });
            foreach my $setup (@setups) {
                my $count = $schema->resultset('DataElementState')->count({ pipelinesetup => $setup->id });
                if ($count < $expected_count) {
                    $rs = $schema->resultset('DataElement')->search({ datasource => $self->id });
                    my @elements;
                    while (my $element = $rs->next) {
                        VRPipe::DataElementState->get(pipelinesetup => $setup, dataelement => $element);
                    }
                }
            }
        }
        
        return 1;
    }
}

1;