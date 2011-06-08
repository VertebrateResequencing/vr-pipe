use VRPipe::Base;

class VRPipe::FileType::lsf with VRPipe::FileTypeRole {
    method _build_read_backwards {
        return 1;
    }
}

1;