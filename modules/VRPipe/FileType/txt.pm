use VRPipe::Base;

class VRPipe::FileType::txt with VRPipe::FileTypeRole {
    method check_type {
        my $file = $self->file;
        return -T $file ? 1 : 0;
    }
}

1;