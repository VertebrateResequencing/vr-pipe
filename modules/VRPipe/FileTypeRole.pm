use VRPipe::Base;

role VRPipe::FileTypeRole {
    use VRPipe::Parser;
    
    has 'file' => (is => 'rw',
                   isa => File,
                   coerce => 1,
                   required => 1);
    
    has 'type' => (is => 'ro',
                   isa => FileType,
                   lazy => 1,
                   builder => '_build_type');
    
    has 'record_separator' => (is => 'ro',
                               isa => 'Maybe[Str]',
                               builder => '_build_record_separator');
    
    has 'read_backwards' => (is => 'ro',
                             isa => 'Bool',
                             builder => '_build_read_backwards');
    
    has 'parser' => (is => 'ro',
                     does => 'VRPipe::ParserRole',
                     lazy => 1,
                     builder => '_build_parser');
    
    method _build_type {
        my $class = ref($self);
        my ($type) = $class =~ /.+::(.+)/;
        return $type;
    }
    
    method _build_record_separator {
        return;
    }
    
    method _build_read_backwards {
        return 0;
    }
    
    method _build_parser {
        return VRPipe::Parser->create($self->type, {file => $self->file});
    }
    
    method check_type {
        my $file = $self->file;
        my $type = $self->type;
        $file =~ s/\.gz$// unless $type eq 'gz';
        if ($file =~ /\.$type$/) {
            return 1;
        }
        return 0;
    }
}

1;