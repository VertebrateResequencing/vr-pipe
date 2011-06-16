use VRPipe::Base;

class VRPipe::DataSource::fofn extends VRPipe::DataSource::list {
    around all (Defined :$handle) {
        my $result;
        
        while ($result = $self->$orig(handle => $handle, skip_comments => 1)) {
            $result = file($result)->absolute; # it would be too expensive to check if this was actually a file by doing -e, so we'll let downstream catch the possible error
            last;
        }
        
        return $result;
    }
}

1;