
=head1 NAME

VRPipe::DataSource::fofn_with_metadata - get pipeline input and meta from a
fofn

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

Similar to L<VRPipe::DataSource::fofn>, the format of the text file is extended
to have an arbitrary number of tab separated columns to hold metadata about
each file, which is then associated with the files in the VRPipe database. The
first line of the file should have the metadata key names. The first column
should always be 'path'.

For example, to use bam files in most pipelines that can process them, the
first line should have the columns: path center_name study sample platform
library lane

This DataSource can also be used to group multiple files together into a single
L<VRPipe::DataElement>, which is essential for 'reduce'/merging type pipelines.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012,2014 Genome Research Limited.

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
        elsif ($method eq 'group_all') {
            return "All files in the file will be grouped into a single element.";
        }
        
        return '';
    }
    
    method all (Defined :$handle!, Str :$filter?) {
        my @element_args;
        my $did = $self->_datasource_id;
        foreach my $result ($self->_all_results(handle => $handle, $filter ? (filter => $filter) : ())) {
            push(@element_args, { datasource => $did, result => { paths => $result->{paths} } });
        }
        $self->_create_elements(\@element_args);
    }
    
    method group_all (Defined :$handle!, Str :$filter?) {
        my @paths;
        foreach my $result ($self->_all_results(handle => $handle, $filter ? (filter => $filter) : ())) {
            push @paths, @{ $result->{paths} };
        }
        return unless @paths;
        $self->_create_elements([{ datasource => $self->_datasource_id, result => { paths => \@paths }, withdrawn => 0 }]);
    }
    
    around _all_results (Defined :$handle!, Str :$delimiter?, ArrayRef :$path_columns?, Bool :$columns_are_paths?, Str :$filter?, Bool :$filter_after_grouping = 0) {
        my @krs;
        if ($filter) {
            foreach my $kr (split(',', $filter)) {
                my ($key, $regex) = split('#', $kr);
                $self->throw("Option 'filter' for vrpipe datasource was not properly formed\n") unless ($key && $regex);
                push(@krs, [$key, $regex]);
            }
            $self->throw("Option 'filter' for vrpipe datasource was not properly formed\n") unless @krs;
        }
        
        my @results;
        my %col_to_key;
        FILE: foreach my $result ($self->$orig(handle => $handle, delimiter => "\t", path_columns => [1])) {
            unless (keys %col_to_key) {
                delete $result->{paths};
                foreach my $col (sort { $a <=> $b } keys %$result) {
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
            my $file = VRPipe::File->create(path => $path);
            $file->add_metadata($metadata, replace_data => 1);
            
            my $pass_filter = 0;
            my $meta        = $file->metadata;
            if ($filter) {
                # if "filter_after_grouping => 0", we filter before grouping
                # by skipping files which don't match the regex or don't
                # have the required metadata
                my $passes = 0;
                foreach my $kr (@krs) {
                    my ($key, $regex) = @$kr;
                    if (defined $meta->{$key}) {
                        my $this_passed = $meta->{$key} =~ m/$regex/ ? 1 : 0;
                        next FILE if (!$filter_after_grouping && !$this_passed);
                        $passes += $this_passed;
                    }
                    else {
                        next FILE unless $filter_after_grouping;
                    }
                }
                $pass_filter = $passes == @krs ? 1 : 0;
            }
            
            push(@results, { paths => [$path], metadata => $metadata });
        }
        
        return @results;
    }
    
    method grouped_by_metadata (Defined :$handle!, Str :$metadata_keys!, Str :$filter?) {
        my @meta_keys = split /\|/, $metadata_keys;
        
        my $group_hash;
        foreach my $result ($self->_all_results(handle => $handle, $filter ? (filter => $filter) : ())) {
            my @group_keys;
            foreach my $key (@meta_keys) {
                $self->throw("Metadata key $key not present for file " . $result->{paths}->[0]) unless (exists $result->{metadata}->{$key});
                push @group_keys, $result->{metadata}->{$key};
            }
            
            my $group_key = join '|', @group_keys;
            push(@{ $group_hash->{$group_key}->{paths} }, @{ $result->{paths} });
        }
        
        my $did = $self->_datasource_id;
        my @element_args;
        while (my ($group, $data) = each %$group_hash) {
            push(@element_args, { datasource => $did, result => { paths => $data->{paths}, group => $group } });
        }
        $self->_create_elements(\@element_args);
    }
}

1;
