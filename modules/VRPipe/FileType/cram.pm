use VRPipe::Base;

class VRPipe::FileType::cram extends VRPipe::FileType::bin {
    our $correct_magic = [qw(103 122 101 115 000 000 000 103 037 213 010 000 000 000 000 000)];
    
    around check_type {
        $self->$orig || return 0;
        my $path = $self->file;
        return $path =~ /\.cram$/ ? 1 : 0;
    }
    
    method num_header_lines (ClassName|Object $self: Str|File :$reference_fasta!) {
        my $path = $self->file;
        my $headers = `java -Dreference=$reference_fasta -cp $ENV{CRAMTOOLS}/cramtools.jar net.sf.picard.sam.ViewSam INPUT=$path | samtools view -SH -`;
        my @header_lines = split(/\n/, $headers);
        return scalar(@header_lines);
    }
    
    method num_records (ClassName|Object $self: Str|File :$reference_fasta!) {
        my $path = $self->file;
        my $records = `java -Dreference=$reference_fasta -cp $ENV{CRAMTOOLS}/cramtools.jar net.sf.picard.sam.ViewSam INPUT=$path | samtools view -Sc -`;
        ($records) = $records =~ /^(\d+)/m;
        $records ||= 0;
        return $records;
    }
}

1;