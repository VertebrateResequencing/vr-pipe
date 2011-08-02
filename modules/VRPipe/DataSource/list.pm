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
    
    method all (Defined :$handle, Bool :$skip_comments = 1, Bool :$line_is_path = 0) {
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
        $result || return;
        
        my $key_name = $line_is_path ? 'paths' : 'line';
        $result = $line_is_path ? [file($result)->absolute->stringify] : $result;
        my %result = ($key_name => $result);
        return \%result;
    }
}

1;