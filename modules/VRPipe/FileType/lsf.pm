use VRPipe::Base;

class VRPipe::FileType::lsf extends VRPipe::FileType::txt {
    method _build_read_backwards {
        return 1;
    }
}

1;