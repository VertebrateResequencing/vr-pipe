use VRPipe::Base;
    
class VRPipe::DataSource::delimited extends VRPipe::DataSource::list {
    method all (Defined :$handle, Str :$delimiter, ArrayRef :$path_columns?) {
        my @elements;
        foreach my $result ($self->_all_results(handle => $handle, delimiter => $delimiter, $path_columns ? ( path_columns => $path_columns) : ())) {
            push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => $result));
        }
        return \@elements;
    }
    
    around _all_results (Defined :$handle, Str :$delimiter, ArrayRef :$path_columns?) {
        my %result;
        $path_columns ||= [];
        my %path_cols = map { $_ => 1 } @$path_columns;
        
        my @results;
        foreach my $result ($self->$orig(handle => $handle, skip_comments => 1)) {
            my @split = split($delimiter, $result->{line});
            
            my $del_result;
            for my $key (1..@split) {
                if (exists $path_cols{$key}) {
                    push(@{$del_result->{paths}}, file($split[$key - 1])->absolute->stringify);
                }
                else {
                    $del_result->{$key} = $split[$key - 1];
                }
            }
            push(@results, $del_result);
        }
        
        return @results;
    }
    
    method single_column (Defined :$handle, Str :$delimiter, PositiveInt :$column, Bool :$column_is_path = 1) {
        my $key_name = $column_is_path ? 'paths' : $column;
        
        my @elements;
        foreach my $hash_ref ($self->_all_results(handle => $handle, delimiter => $delimiter, $column_is_path ? (path_columns => [$column]) : ())) {
            push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => { $key_name => $hash_ref->{$key_name} }));
        }
        
        return \@elements;
    }
    
    method grouped_single_column (Defined :$handle, Str :$delimiter, PositiveInt :$column, PositiveInt :$group_by, Bool :$column_is_path = 1) {
        my $group_hash;
        foreach my $hash_ref ($self->_all_results(handle => $handle, delimiter => $delimiter, $column_is_path ? (path_columns => [$column]) : ())) {
            push(@{$group_hash->{$hash_ref->{$group_by}}}, $column_is_path ? @{$hash_ref->{paths}} : $hash_ref->{$column});
        }
        
        my $key_name = $column_is_path ? 'paths' : $column;
        my @elements;
        foreach my $group (sort keys %$group_hash) {
            my $array_ref = $group_hash->{$group};
            push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => { $key_name => $array_ref, group => $group }, withdrawn => 0));
        }
        
        return \@elements;
    }
}

1;