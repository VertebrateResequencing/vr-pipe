use VRPipe::Base;

class VRPipe::FileType::bcf extends VRPipe::FileType::bin {
    our $correct_magic     = [qw(102 103 106 002)];
    our $old_correct_magic = [qw(037 213 010 004 000 000 000 000 000 377 006 000 102 103 002)];
    
    around check_type {
        $self->$orig || return 0;
        my $file = $self->file;
        $self->check_magic($file, $correct_magic) || $self->check_magic($file, $old_correct_magic) || return 0;
        
        # unfortunately bgzip and bams have the same magic, so we'll incorrectly
        # validate compressed VCFs. One way to get around this is to decompress
        # a little bit and see if the BCF\1 magic is there, but we'll just
        # enforce that the filename ends in .bcf for now
        if ($file =~ /\.bcf$/) {
            return 1;
        }
        return 0;
    }
    
    method num_header_lines {
        return scalar(@{ $self->_header_lines });
    }
    
    method num_records {
        my $path    = $self->file;
        my $lines   = `bcftools view $path | wc -l`;
        my $records = $lines - $self->num_header_lines();
        return $records;
    }
    
    method samples {
        my @line = split(/\t/, ${ $self->_header_lines }[-1]);
        my @samples = @line[9 .. $#line];
        @samples || $self->throw("No samples found in bcf, " . $self->file);
        return \@samples;
    }
    
    method _header_lines {
        my @header_lines;
        my $bcf  = $self->file;
        my $open = "bcftools view $bcf |";
        open(my $fh, $open) || $self->throw("Couldn't open '$open': $!");
        while (<$fh>) {
            chomp;
            if (/^#/) {
                push(@header_lines, $_);
            }
            else {
                $self->throw("No header line found for $bcf");
            }
            last if /^#CHROM/;
        }
        close($fh);
        return \@header_lines;
    }
}

1;
