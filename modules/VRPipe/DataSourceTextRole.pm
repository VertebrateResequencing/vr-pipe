use VRPipe::Base;

role VRPipe::DataSourceTextRole with VRPipe::DataSourceRole {
    has 'source_file' => => (is => 'ro',
                             isa => 'VRPipe::File',
                             lazy => 1,
                             builder => '_build_source_file');
    
    method _build_source_file {
        my $source = file($self->source)->absolute;
        return VRPipe::File->get(path => $source, type => 'txt');
    }
    
    method _has_changed {
        my $old_md5 = $self->_changed_marker || return 1;
        my $file = $self->source_file;
        my $new_md5 = $self->file_md5($file);
        if ($new_md5 ne $old_md5) {
            $file->update_md5($new_md5);
            return 1;
        }
        else {
            return 0;
        }
    }
    
    method _update_changed_marker {
        my $file = $self->source_file;
        my $md5 = $self->file_md5($file);
        $file->update_md5($md5);
        $self->_changed_marker($md5);
    }
}

1;