use VRPipe::Base;

class VRPipe::FileType::cat extends VRPipe::FileType::txt {
    method _build_read_backwards {
        return 1;
    }
}

1;