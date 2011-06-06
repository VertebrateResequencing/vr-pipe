use VRPipe::Base;

class VRPipe::FileDefinition extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'type' => (is => 'rw',
                   isa => FileType,
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'match_sub' => (is => 'rw',
                        isa => 'CodeRef',
                        traits => ['VRPipe::Persistent::Attributes'],
                        builder => '_default_match_sub');
    
    has 'output_sub' => (is => 'rw',
                         isa => 'CodeRef',
                         traits => ['VRPipe::Persistent::Attributes'],
                         builder => '_default_output_sub');
    
    has 'description' => (is => 'rw',
                          isa => Varchar[64],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1);
    
    method _default_match_sub {
        sub {
            my ($self, $file) = @_;
            my ($type) = $file =~ /\.(\w{1,3})$/;
            if ($type) {
                return $type eq $self->type ? 1 : 0;
            }
            return 0;
        }
    }
    
    method _default_output_sub {
        sub {
            my $self = shift;
            my $basename = $self->name.'.output.'.$self->type;
            $basename =~ s/\W/_/;
            return $basename;
        }
    }
    
    __PACKAGE__->make_persistent();
    
    method matches (VRPipe::File|Str|File $file) {
        # convert input to a Path::Class::File
        if (ref($file)) {
            if ($file->isa('VRPipe::File')) {
                return 0 unless $self->type eq $file->type;
                $file = $file->path;
            }
        }
        else {
            $file = file($file);
        }
        
        my $sub = $self->match_sub();
        return &$sub($self, $file);
    }
    
    method output_basename (PersistentFileHashRef :$inputs, HashRef :$options, VRPipe::DataElement :$data_element) {
        my $sub = $self->output_sub();
        return &$sub($self, $inputs, $options, $data_element);
    }
}

1;