use VRPipe::Base;

class VRPipe::FileType::bcf extends VRPipe::FileType::vcf {
    method check_type {
        return ($self->hts_file_type =~ /^BCF/) ? 1 : 0;
    }

}

1;
