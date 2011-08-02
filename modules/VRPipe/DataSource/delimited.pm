use VRPipe::Base;
    
class VRPipe::DataSource::delimited extends VRPipe::DataSource::list {
    has '_group_hash' => (is => 'rw',
                          isa => 'HashRef');
    
    around all (Defined :$handle, Str :$delimiter, ArrayRef :$path_columns?) {
        my %result;
        $path_columns ||= [];
        my %path_cols = map { $_ => 1 } @$path_columns;
        
        while (my $result = $self->$orig(handle => $handle, skip_comments => 1)) {
            my @split = split($delimiter, $result->{line});
            
            for my $key (1..@split) {
                if (exists $path_cols{$key}) {
                    push(@{$result{paths}}, file($split[$key - 1])->absolute->stringify);
                }
                else {
                    $result{$key} = $split[$key - 1];
                }
            }
            
            last;
        }
        
        return keys %result ? \%result : undef;
    }
    
    method single_column (Defined :$handle, Str :$delimiter, PositiveInt :$column, Bool :$column_is_path = 1) {
        my %result;
        
        while (my $hash_ref = $self->all(handle => $handle, delimiter => $delimiter, $column_is_path ? (path_columns => [$column]) : ())) {
            my $key_name = $column_is_path ? 'paths' : $column;
            $result{$key_name} = $hash_ref->{$key_name};
            last;
        }
        
        return keys %result ? \%result : undef;
    }
    
    method grouped_single_column (Defined :$handle, Str :$delimiter, PositiveInt :$column, PositiveInt :$group_by, Bool :$column_is_path = 1) {
        my $group_hash = $self->_group_hash;
        
        unless ($group_hash) {
            while (my $hash_ref = $self->all(handle => $handle, delimiter => $delimiter, $column_is_path ? (path_columns => [$column]) : ())) {
                push(@{$group_hash->{$hash_ref->{$group_by}}}, $column_is_path ? @{$hash_ref->{paths}} : $hash_ref->{$column});
            }
            $self->_group_hash($group_hash);
        }
        
        my @groups = sort keys %$group_hash;
        @groups || return;
        
        my $group = $groups[0];
        my $array_ref = delete $group_hash->{$group};
        my $key_name = $column_is_path ? 'paths' : $column;
        my %result = ($key_name => $array_ref,
                      group => $group);
        return \%result;
    }
}

1;