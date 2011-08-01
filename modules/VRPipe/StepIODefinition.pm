use VRPipe::Base;

class VRPipe::StepIODefinition extends VRPipe::Persistent {
    has 'type' => (is => 'rw',
                   isa => FileType,
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'min_files' => (is => 'rw',
                        isa => IntSQL[4],
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_key => 1,
                        default => 1,
                        allow_key_to_default => 1);
    
    has 'max_files' => (is => 'rw',
                        isa => IntSQL[4],
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_key => 1,
                        default => 1, # -1 means no maximum
                        allow_key_to_default => 1);
    
    has 'metadata' => (is => 'rw',
                       isa => 'HashRef',
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       default => sub { {} },
                       allow_key_to_default => 1);
    
    has 'description' => (is => 'rw',
                          isa => Varchar[255],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1);
    
    __PACKAGE__->make_persistent();
    
    method required_metadata_keys {
        my $metadata = $self->metadata;
        my %optional = map { $_ => 1 } @{delete $metadata->{optional} || []};
        
        my @required;
        foreach my $key (sort keys %$metadata) {
            next if exists $optional{$key};
            push(@required, $key);
        }
        
        return @required;
    }
}

1;