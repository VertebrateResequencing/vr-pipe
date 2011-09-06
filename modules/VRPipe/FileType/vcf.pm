use VRPipe::Base;

class VRPipe::FileType::vcf extends VRPipe::FileType::txt {
    method check_type {
        # we don't use txt type check, because we allow for compressed and
        # uncompressed vcf files
        
        #*** worth doing something like checking the first line of file?
        
        return 1;
    }
    
    method num_header_lines {
        my $path = $self->file;
        my $vrfile = VRPipe::File->get(path => $path);
        my $fh = $vrfile->openr;
        my $count = 0;
        while (<$fh>) {
            if (/^#/) {
                $count++;
            }
            else {
                last;
            }
        }
        $vrfile->close;
        return $count;
    }
    around num_records {
        my $total_lines = $self->$orig();
        my $header_lines = $self->num_header_lines;
        return $total_lines - $header_lines;
    }
}

1;