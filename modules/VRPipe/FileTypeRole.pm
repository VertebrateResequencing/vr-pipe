use VRPipe::Base;

role VRPipe::FileTypeRole {
    has 'file' => (is => 'rw',
                   isa => File,
                   coerce => 1,
                   required => 1);
    
    has 'type' => (is => 'ro',
                   isa => FileType,
                   lazy => 1,
                   builder => '_build_type');
    
    method _build_type {
        my $class = ref($self);
        my ($type) = $class =~ /.+::(.+)/;
        return $type;
    }
    
    method check_type {
        my $file = $self->file;
        my $type = $self->type;
        if ($file =~ /\.$type$/) {
            return 1;
        }
        return 0;
    }
}

1;