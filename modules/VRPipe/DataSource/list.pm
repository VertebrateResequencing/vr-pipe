use VRPipe::Base;

class VRPipe::DataSource::list with VRPipe::DataSourceRole {
    use VRPipe::File;
    
    has '_vrpfile' => (is => 'rw',
                       isa => 'VRPipe::File');
    
    method _open_source {
        my $source = file($self->source)->absolute;
        my $file = VRPipe::File->get(path => $source, type => 'txt');
        $file->e || $self->throw("list file $source does not exist!");
        $self->_vrpfile($file); # so that it isn't destroyed, which would close the handle
        return $file->openr;
    }
    
    method all (Defined :$handle, Bool :$skip_comments = 1) {
        my $result;
        
        while (<$handle>) {
            if (/^\s*$/) {
                next;
            }
            elsif ($skip_comments && /^#/) {
                next;
            }
            chomp;
            $result = $_;
            last;
        }
        
        return $result;
    }
}

1;