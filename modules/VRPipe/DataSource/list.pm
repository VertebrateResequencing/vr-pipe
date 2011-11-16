use VRPipe::Base;

class VRPipe::DataSource::list with VRPipe::DataSourceTextRole {
    use VRPipe::File;
    
    method description {
        return "Use a simple list of items in a file as your source.";
    }
    method source_description {
        return "The path to a file with one item per line.";
    }
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "Each element will correspond to a single line from the file.";
        }
        
        return '';
    }
    
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