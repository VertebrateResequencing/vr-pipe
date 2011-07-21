use VRPipe::Base;

class VRPipe::FileType::bam extends VRPipe::FileType::bin {
    our $correct_magic = [qw(037 213 010 004 000 000 000 000 000 377 006 000 102 103 002 000)];
    
    around check_type {
        $self->$orig || return 0;
        return $self->check_magic($self->file, $correct_magic);
    }
    
    method num_header_lines {
        my $path = $self->file;
        my $headers = `samtools view -H $path`;
        my @header_lines = split(/\n/, $headers);
        return scalar(@header_lines);
    }
    method num_records {
        my $path = $self->file;
        my $records = `samtools view -c $path`;
        ($records) = $records =~ /^(\d+)/m;
        return $records;
    }
}

1;