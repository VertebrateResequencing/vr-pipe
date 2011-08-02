use VRPipe::Base;

class VRPipe::FileType::any with VRPipe::FileTypeRole {
    method check_type {
        return 1;
    }
}

1;