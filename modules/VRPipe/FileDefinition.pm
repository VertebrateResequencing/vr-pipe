use VRPipe::Base;

class VRPipe::FileDefinition extends VRPipe::Persistent {
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
                          isa => Varchar[256],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1);
    
    method _default_match_sub {
        sub {
            my ($self, $file) = @_;
            my $type = VRPipe::FileType->create($self->type, {file => $file});
            return $type->check_type;
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
                # return 0 unless $self->type eq $file->type; *** this doesn't allow subtypes (eg. 'txt') to match supertypes (eg. 'any'); is there an alternative shortcut we can do?
                $file = $file->path;
            }
        }
        else {
            $file = file($file);
        }
        
        my $sub = $self->match_sub();
        return &$sub($self, $file);
    }
    
    method output_basename (VRPipe::Step $step) {
        my $sub = $self->output_sub();
        return &$sub($self, $step);
    }
}

1;