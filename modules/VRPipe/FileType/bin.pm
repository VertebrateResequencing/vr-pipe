use VRPipe::Base;

class VRPipe::FileType::bin with VRPipe::FileTypeRole {
    method check_type {
        my $file = $self->file;
        return -B $file ? 1 : 0;
    }
}

1;