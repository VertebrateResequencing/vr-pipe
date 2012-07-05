=head1 NAME

VRPipe::DataSource::delimited - get pipeline input from a delimited text file

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This DataSource has a variety of methods to work with tab-delimited text files.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;
    
class VRPipe::DataSource::delimited extends VRPipe::DataSource::list {
    method description {
        return "Use information in a delimited text file.";
    }
    method source_description {
        return "The path to a delimited text file.";
    }
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "Each element will consist of a list of all the fields on a single line.";
        }
        elsif ($method eq 'single_column') {
            return "Each element will consist of the entry of a particular column from a single line.";
        }
        elsif ($method eq 'all_columns') {
            return "Each element will consist of a list of all the file paths on a single line.";
        }
        elsif ($method eq 'grouped_single_column') {
            return "Each element will consist of the file paths given in a particular column, grouped by another column.";
        }
        
        return '';
    }
    
    method all (Defined :$handle!, Str :$delimiter!, ArrayRef :$path_columns?) {
        my @element_args;
        foreach my $result ($self->_all_results(handle => $handle, delimiter => $delimiter, $path_columns ? ( path_columns => $path_columns) : ())) {
            push(@element_args, { datasource => $self->_datasource_id, result => $result });
        }
        $self->_create_elements(\@element_args);
    }
    
    around _all_results (Defined :$handle!, Str :$delimiter!, ArrayRef :$path_columns?, Bool :$columns_are_paths?) {
        my %result;
        $path_columns ||= [];
        my %path_cols = map { $_ => 1 } @$path_columns;
        
        my @results;
        foreach my $result ($self->$orig(handle => $handle, skip_comments => 1)) {
            my @split = split($delimiter, $result->{line});
            
            my $del_result;
            for my $key (1..@split) {
                if ($columns_are_paths || exists $path_cols{$key}) {
                    push(@{$del_result->{paths}}, VRPipe::File->create(path => file($split[$key - 1])->absolute)->path->stringify); # we can't bulk_create VRPipe::Files because they do fancy stuff duing create()
                }
                else {
                    $del_result->{$key} = $split[$key - 1];
                }
            }
            push(@results, $del_result);
        }
        
        return @results;
    }
    
    method single_column (Defined :$handle!, Str :$delimiter!, PositiveInt :$column!, Bool :$column_is_path = 1) {
        my $key_name = $column_is_path ? 'paths' : $column;
        
        my @element_args;
        foreach my $hash_ref ($self->_all_results(handle => $handle, delimiter => $delimiter, $column_is_path ? (path_columns => [$column]) : ())) {
            push(@element_args, { datasource => $self->_datasource_id, result => { $key_name => $hash_ref->{$key_name} } });
        }
        $self->_create_elements(\@element_args);
    }
    
    method all_columns (Defined :$handle!, Str :$delimiter!) {
        my @element_args;
        foreach my $hash_ref ($self->_all_results(handle => $handle, delimiter => $delimiter, columns_are_paths => 1)) {
            push(@element_args, { datasource => $self->_datasource_id, result => { paths => $hash_ref->{paths} }});
        }
        $self->_create_elements(\@element_args);
    }
    
    method grouped_single_column (Defined :$handle!, Str :$delimiter!, PositiveInt :$column!, PositiveInt :$group_by!, Bool :$column_is_path = 1) {
        my $group_hash;
        foreach my $hash_ref ($self->_all_results(handle => $handle, delimiter => $delimiter, $column_is_path ? (path_columns => [$column]) : ())) {
            push(@{$group_hash->{$hash_ref->{$group_by}}}, $column_is_path ? @{$hash_ref->{paths}} : $hash_ref->{$column});
        }
        
        my $key_name = $column_is_path ? 'paths' : $column;
        my @element_args;
        foreach my $group (sort keys %$group_hash) {
            my $array_ref = $group_hash->{$group};
            push(@element_args, { datasource => $self->_datasource_id, result => { $key_name => $array_ref, group => $group } });
        }
        $self->_create_elements(\@element_args);
    }
}

1;