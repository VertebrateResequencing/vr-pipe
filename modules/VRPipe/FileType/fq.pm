use VRPipe::Base;

class VRPipe::FileType::fq extends VRPipe::FileType::txt {
    #*** needs more implementation!
    
    method check_type {
        my $file = $self->file;
        my $type = $self->type;
        $file =~ s/\.gz$// unless $type eq 'gz';
        if ($file =~ /\.(?:$type|fastq)$/) { #*** this sucks as a test...
            return 1;
        }
        return 0;
    }
}

1;