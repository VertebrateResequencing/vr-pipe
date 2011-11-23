use VRPipe::Base;

class VRPipe::DataSource::fofn extends VRPipe::DataSource::list {
    method description {
        return "Use a simple list of absolute file paths (file-of-file-names) in a file as your source.";
    }
    method source_description {
        return "The path to a file with one absolute file path per line.";
    }
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "Each element will correspond to a single file from the file.";
        }
        
        return '';
    }
    
    method all (Defined :$handle!) {
        my @elements;
        foreach my $result ($self->_all_results(handle => $handle, skip_comments => 1, line_is_path => 1)) {
            my $withdraw = -e $result->{paths}->[0] ? 0 : 1;
            push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => $result, withdrawn => $withdraw));
        }
        return \@elements;
    }
}

1;