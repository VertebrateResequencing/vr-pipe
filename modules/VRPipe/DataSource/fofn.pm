use VRPipe::Base;

class VRPipe::DataSource::fofn extends VRPipe::DataSource::list {
    around all (Defined :$handle) {
        # it would be too expensive to check if each line of the file was
        # actually a file by doing -e, so we'll let downstream catch the
        # possible error
        
        return $self->$orig(handle => $handle, skip_comments => 1, line_is_path => 1);
    }
}

1;