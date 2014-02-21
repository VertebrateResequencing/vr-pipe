
=head1 NAME

VRPipe::DataSource::irods - get pipeline inputs from an iRods data store

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013,2014 Genome Research Limited.

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

class VRPipe::DataSource::irods with VRPipe::DataSourceRole {
    use Digest::MD5 qw(md5_hex);
    use File::Spec::Functions;
    use VRPipe::Persistent::InMemory;
    
    has '_irods_files_and_metadata_cache' => (
        is        => 'rw',
        isa       => 'HashRef',
        clearer   => '_clear_cache',
        predicate => '_cached',
    );
    
    method description {
        return "Use an iRods data store to extract information from";
    }
    
    method source_description {
        return "The zone to query";
    }
    
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "An element will comprise one of the files returned by imeta qu -d given the arguments you supply for the 'file_query' option. The file will have all the relevant irods metadata associated with it, and a local path based on the 'local_root_dir' option. To avoid spamming the irods server, the update_interval option allows you to specify the minimum number of minutes between each check for changes to files. If update_interval is not supplied it defaults to 1 minute when testing, and 1 day in production. The Sanger-specific (otherwise ignorable) add_metadata_from_warehouse option (enabled by default, requiring the environment variables WAREHOUSE_DATABASE, WAREHOUSE_HOST, WAREHOUSE_PORT and WAREHOUSE_USER) will add extra metadata to the files with the key public_name (if defined).";
        }
        return '';
    }
    
    method _open_source {
        return $self->source; # the source can't really be opened
    }
    
    method _has_changed {
        my $old      = $self->_changed_marker;
        my $checksum = $self->_irods_file_metadata_checksum;
        $self->_changed_marker($checksum);
        $old || return 1;
        return ($checksum ne $old) ? 1 : 0;
    }
    
    method _update_changed_marker {
        $self->_changed_marker($self->_irods_file_metadata_checksum);
    }
    
    method _irods_file_metadata_checksum {
        # we need a quick way of looking at all the relevant metadata for all
        # our desired files, but imeta/iquest doesn't give us this; to avoid
        # abusing the irods server too much we only actually calculate the
        # checksum anew every x mins; otherwise we just return the current
        # checksum
        my $im              = VRPipe::Persistent::InMemory->new;
        my $options         = $self->options;
        my $update_interval = $options->{update_interval} || $VRPipe::Persistent::InMemory::deployment eq 'production' ? 1440 : 0; # in mins
        $update_interval *= 60;                                                                                                    # in seconds
        $update_interval ||= 5;                                                                                                    # for testing
        my $lock_key         = 'irods_datasource.' . $self->_datasource_id;
        my $locked           = $im->lock($lock_key, unlock_after => $update_interval);
        my $current_checksum = $self->_changed_marker;
        
        if ($current_checksum && length($current_checksum) == 32) {
            return $current_checksum unless $locked;
        }
        # else we always get the latest checksum if we have no valid checksum
        
        # get the current files and their metadata and stringify it all
        my $files = $self->_get_irods_files_and_metadata;
        $self->_irods_files_and_metadata_cache($files);
        my $data = '';
        foreach my $file (sort keys %$files) {
            $data .= $file;
            
            my $meta = $files->{$file};
            foreach my $key (sort keys %$meta) {
                next if $key eq 'vrpipe_irods_order';
                my $val = $meta->{$key};
                $data .= "|$key:$val|";
            }
        }
        
        my $digest = md5_hex $data;
        return $digest;
    }
    
    method _get_irods_files_and_metadata (Str $zone?, Str $query?) {
        return $self->_irods_files_and_metadata_cache if $self->_cached;
        $zone  ||= $self->_open_source();
        $query ||= $self->options->{file_query};
        
        my $cmd = "imeta -z $zone qu -d $query";
        open(my $qu_fh, $cmd . ' |') || $self->throw("Could not open pipe to [$cmd]");
        my %files;
        my $collection;
        my $order = 1;
        while (<$qu_fh>) {
            #*** do we have to worry about spaces in file paths?...
            if (/^collection:\s+(\S+)/) {
                $collection = $1;
            }
            elsif (/^dataObj:\s+(\S+)/) {
                my $path = "$collection/$1";
                
                # get all the metadata for this file
                my $ls_cmd = "imeta ls -d $path";
                open(my $ls_fh, $ls_cmd . ' |') || $self->throw("Could not open pipe to [$ls_cmd]");
                my $meta = { vrpipe_irods_order => $order++ };
                my $attribute;
                while (<$ls_fh>) {
                    if (/^attribute:\s+(\S+)/) {
                        $attribute = $1;
                        undef $attribute if $attribute =~ /^dcterms:/;
                    }
                    elsif ($attribute && /^value:\s+(.+)$/) {
                        $meta->{$attribute} = $1;
                    }
                }
                close($ls_fh);
                
                $files{$path} = $meta;
            }
        }
        close($qu_fh);
        
        return \%files;
    }
    
    method all (Defined :$handle!, Str :$file_query!, Str|Dir :$local_root_dir!, Str :$update_interval?, Bool :$add_metadata_from_warehouse = 1) {
        # _get_irods_files_and_metadata will get called twice in row: once to
        # see if the datasource changed, and again here; _has_changed caches
        # the result, and we clear the cache after getting that data
        my $files = $self->_get_irods_files_and_metadata($handle, $file_query);
        $self->_clear_cache;
        
        my ($warehouse_sth, $public_name);
        if ($add_metadata_from_warehouse && $ENV{WAREHOUSE_DATABASE} && $ENV{WAREHOUSE_HOST} && $ENV{WAREHOUSE_PORT} && $ENV{WAREHOUSE_USER}) {
            my $dbh = DBI->connect(
                "DBI:mysql:host=$ENV{WAREHOUSE_HOST}:port=$ENV{WAREHOUSE_PORT};database=$ENV{WAREHOUSE_DATABASE}",
                $ENV{WAREHOUSE_USER}, undef, { 'RaiseError' => 1, 'PrintError' => 0 }
            );
            
            # sanger_sample_id in warehouse corresponds to 'sample' metadata
            # from irods, and is one of the few indexed columns, so queries
            # against it should hopefully be quick. (supplier_name in warehouse
            # is sample_supplier_name in irods, and donor_id in warehouse is
            # sample_cohort in irods, we don't need to also get those columns)
            my $sql = q[select public_name from current_samples where sanger_sample_id = ?];
            
            $warehouse_sth = $dbh->prepare($sql);
            $warehouse_sth->execute;
            $warehouse_sth->bind_col(1, \$public_name);
        }
        
        my $did = $self->_datasource_id;
        my @element_args;
        foreach my $path (sort { $files->{$a}->{vrpipe_irods_order} <=> $files->{$b}->{vrpipe_irods_order} } keys %$files) {
            my $new_metadata = $files->{$path};
            delete $new_metadata->{vrpipe_irods_order};
            if ($warehouse_sth) {
                undef $public_name;
                $warehouse_sth->execute($new_metadata->{sample});
                $warehouse_sth->fetch;
                if ($public_name) {
                    $new_metadata->{public_name} = "$public_name";
                }
            }
            
            my $sub_path = $path;
            $sub_path =~ s/^\///;
            my $file_abs_path = file($local_root_dir, $sub_path)->stringify;
            
            # consider type to be txt if not a VRPipe file type
            my $type = delete $new_metadata->{type};
            $type ||= 'txt';
            eval "require VRPipe::FileType::$type;";
            if ($@) { $type = 'txt'; }
            
            my $vrfile = VRPipe::File->create(path => $file_abs_path, type => $type)->original;
            
            # add metadata to file, detecting any changes
            my $current_metadata = $vrfile->metadata;
            my $changed          = 0;
            if ($current_metadata && keys %$current_metadata) {
                while (my ($key, $val) = each %$current_metadata) {
                    next unless defined $val;
                    next unless defined $new_metadata->{$key};
                    if ($val ne $new_metadata->{$key}) {
                        $self->debug("metadata '$key' changed from $val to $new_metadata->{$key} for file $file_abs_path, so will mark file as changed");
                        $changed = 1;
                        last;
                    }
                }
            }
            
            # if there was no metadata this will add metadata to the file.
            $vrfile->add_metadata($new_metadata, replace_data => 0);
            $vrfile->add_metadata({ irods_path => $path });
            
            my $result_hash = { paths => [$file_abs_path], irods_path => $path };
            if ($changed) {
                $result_hash->{changed} = [[$vrfile, $new_metadata]];
                $self->_start_over_elements_due_to_file_metadata_change($result_hash);
                delete $result_hash->{changed};
            }
            push(@element_args, { datasource => $did, result => $result_hash });
        }
        $self->_create_elements(\@element_args);
    }
}

1;
