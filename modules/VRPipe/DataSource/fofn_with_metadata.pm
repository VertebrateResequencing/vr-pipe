use VRPipe::Base;
    
class VRPipe::DataSource::fofn_with_metadata extends VRPipe::DataSource::delimited {
    method description {
        return "Use information in a tab-delimited text file which specifies file paths and metadata to apply to that file.";
    }
    method source_description {
        return "The path to a tab-delimeted text file where the first line has columns 'path' and then metadata key names, and subsequent lines have the absolute path to a file and the corresponding metadata key values. For example, to use bam files in most pipelines, the first line should have the columns: path center_name study sample platform library lane.";
    }
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "Each element will consist of a single file given in the first column of the source; subsequent columns are applied as metadata to the file (overriding any existing metadata).";
        }
        elsif ($method eq 'grouped_by_metadata') {
            return "Each element will consist of all the files given in the first column of the source that share values in the column(s) specified by the metadata_keys option. If you want to group by more than 1 metadata column, separate the key names with a '|' symbol. e.g. 'sample|platform|library' will group all files with the same sample, platform and library into one dataelement.";
        }
        
        return '';
    }
    
    method all (Defined :$handle!) {
        my @elements;
        foreach my $result ($self->_all_results(handle => $handle)) {
            push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => {paths => $result->{paths}}));
        }
        return \@elements;
    }
    
    around _all_results (Defined :$handle!, Str :$delimiter?, ArrayRef :$path_columns?, Bool :$columns_are_paths?) {
        my @results;
        my %col_to_key;
        foreach my $result ($self->$orig(handle => $handle, delimiter => "\t", path_columns => [1])) {
            unless (keys %col_to_key) {
                delete $result->{paths};
                foreach my $col (sort {$a <=> $b} keys %$result) {
                    $col_to_key{$col} = $result->{$col};
                }
                unless (keys %col_to_key) {
                    $self->throw("source had no metadata columns?");
                }
                next;
            }
            
            my $metadata = {};
            while (my ($col, $key) = each %col_to_key) {
                $metadata->{$key} = $result->{$col};
            }
            
            my $path = file($result->{paths}->[0])->absolute->stringify;
            my $file = VRPipe::File->get(path => $path);
            $file->add_metadata($metadata, replace_data => 1);
            
            push(@results, {paths => [$path], metadata => $metadata});
        }
        
        return @results;
    }
    
    method grouped_by_metadata (Defined :$handle!, Str :$metadata_keys!) {
        my @meta_keys = split /\|/, $metadata_keys;
        
        my $group_hash;
        foreach my $result ($self->_all_results(handle => $handle)) {
            my @group_keys;
            foreach my $key (@meta_keys) {
                $self->throw("Metadata key $key not present for file ".$result->{paths}->[0]) unless (exists $result->{metadata}->{$key});
                push @group_keys, $result->{metadata}->{$key};
            }
            
            my $group_key = join '|', @group_keys;
            push(@{$group_hash->{$group_key}->{paths}}, @{$result->{paths}});
        }
        
        my $did = $self->_datasource_id;
        my @elements;
        while (my ($group, $data) = each %$group_hash) {
            push(@elements, VRPipe::DataElement->get(datasource => $did, result => { paths => $data->{paths}, group => $group }, withdrawn => 0));
        }
        return \@elements;
    }
}

1;