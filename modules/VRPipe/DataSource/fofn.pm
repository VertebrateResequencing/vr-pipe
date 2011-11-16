use VRPipe::Base;

class VRPipe::DataSource::fofn extends VRPipe::DataSource::list {
    method description {
        return "Use a simple list of absolute file paths (file-of-file-names) in a file as your source.";
    }
    method source_description {
        return "The path to a file with one absolute file path per line.";
    }
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "Each element will correspond to a single file from the file.";
        }
        
        return '';
    }
    
    around all (Defined :$handle) {
        # it would be too expensive to check if each line of the file was
        # actually a file by doing -e, so we'll let downstream catch the
        # possible error
        
        return $self->$orig(handle => $handle, skip_comments => 1, line_is_path => 1);
    }
}

1;