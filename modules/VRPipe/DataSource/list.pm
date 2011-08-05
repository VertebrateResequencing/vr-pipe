use VRPipe::Base;

class VRPipe::DataSource::list with VRPipe::DataSourceTextRole {
    use VRPipe::File;
    
    method _open_source {
        my $file = $self->source_file;
        return $file->openr;
    }
    
    method all (Defined :$handle, Bool :$skip_comments = 1, Bool :$line_is_path = 0) {
        my @elements;
        foreach my $result ($self->_all_results(handle => $handle, skip_comments => $skip_comments, line_is_path => $line_is_path)) {
            push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => $result, withdrawn => 0));
        }
        return \@elements;
    }
    
    method _all_results (Defined :$handle, Bool :$skip_comments = 1, Bool :$line_is_path = 0) {
        my $key_name = $line_is_path ? 'paths' : 'line';
        
        my @results;
        while (<$handle>) {
            if (/^\s*$/) {
                next;
            }
            elsif ($skip_comments && /^#/) {
                next;
            }
            chomp;
            
            my $result = $_;
            $result = $line_is_path ? [file($result)->absolute->stringify] : $result;
            push(@results, { $key_name => $result });
        }
        
        return @results;
    }
}

1;