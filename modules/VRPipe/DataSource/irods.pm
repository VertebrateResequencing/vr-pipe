
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
    use Path::Class;
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
            return "An element will comprise one of the files returned by imeta qu -d given the arguments you supply for the 'file_query' option. The file will have all the relevant irods metadata associated with it, and a local path based on the 'local_root_dir' option. To avoid spamming the irods server, the update_interval option allows you to specify the minimum number of minutes between each check for changes to files. If update_interval is not supplied it defaults to 5 seconds when testing, and 1 day in production.";
        }
        if ($method eq 'all_with_warehouse_metadata') {
            return "In addition to doing everything the all method does, it adds extra metadata found in the warehouse database to the files with the keys public_name, sample_supplier_name, sample_control, sample_cohort, taxon_id, sample_created_date and study_title (if defined). If any analysis has been done to a file, the associated files are stored under irods_analysis_files. (This method is Sanger-specific and also requires the environment variables WAREHOUSE_DATABASE, WAREHOUSE_HOST, WAREHOUSE_PORT and WAREHOUSE_USER.)";
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
        my $files = $self->_get_irods_files_and_metadata($self->_open_source(), $options->{file_query}, $self->method eq 'all_with_warehouse_metadata');
        $self->_irods_files_and_metadata_cache($files);
        my $data = '';
        foreach my $file (sort keys %$files) {
            $data .= $file;
            
            my $meta = $files->{$file};
            foreach my $key (sort keys %$meta) {
                next if $key eq 'vrpipe_irods_order';
                my $val = $meta->{$key};
                # limit to printable ascii so md5_hex subroutine will work
                $val =~ tr/\x20-\x7f//cd;
                $data .= "|$key:$val|";
            }
        }
        
        my $digest = md5_hex $data;
        return $digest;
    }
    
    method _get_irods_files_and_metadata (Str $zone!, Str $query!, Bool $add_metadata_from_warehouse?) {
        return $self->_irods_files_and_metadata_cache if $self->_cached;
        
        my ($sample_sth, $public_name, $donor_id, $supplier_name, $control, $taxon_id, $created);
        my ($study_sth, $study_title);
        if ($add_metadata_from_warehouse && $ENV{WAREHOUSE_DATABASE} && $ENV{WAREHOUSE_HOST} && $ENV{WAREHOUSE_PORT} && $ENV{WAREHOUSE_USER}) {
            my $dbh = DBI->connect(
                "DBI:mysql:host=$ENV{WAREHOUSE_HOST}:port=$ENV{WAREHOUSE_PORT};database=$ENV{WAREHOUSE_DATABASE}",
                $ENV{WAREHOUSE_USER}, undef, { 'RaiseError' => 1, 'PrintError' => 0, mysql_enable_utf8 => 1 }
            );
            
            # name in warehouse corresponds to 'sample' metadata from irods, and
            # is one of the few indexed columns, so queries against it should
            # hopefully be quick. For situations where the irods metadata is
            # missing or out-of-date we want to get values from warehouse and
            # store them under the metadata key that irods uses:
            # warehouse_table,warehouse_column => irods_key
            # current_samples,supplier_name => sample_supplier_name
            # current_samples,donor_id => sample_cohort
            # current_samples,control => sample_control
            # current_samples,taxon_id => n/a (store under taxon_id)
            # current_samples,created => n/a (store under sample_created_date)
            # current_studies,name => study_title (linked to irods study_id => internal_id)
            
            # notes for /archive/GAPI/exp/infinium idat files:
            # lane name = basename of filename
            # library name = sample.beadchip.beadchip_section
            # library ssid = last 4 digits of beadchip . last 3 digits of sample_id
            # lane acc = analysis_uuid (the last one listed by imeta)
            # lane storage_path = {user dir}/analysis_uuid
            
            # notes for /archive/GAPI/gen/infinium gtc files:
            # lane name = basename of filename, or {beadchip}_{beadchip_section}
            # library name = infinium_sample
            # library ssid = last 4 digits of beadchip . last 3 digits of sample_id
            # lane acc = analysis_uuid (the last one listed by imeta)
            # lane storage_path = {user dir}/{analysis_uuid}_{basename(user dir)}
            
            # notes for /seq bam files:
            # lane name = basename of filename
            # individual name = sample_supplier_name or sample
            # individual acc = sample_accession_number
            
            my $sql = q[select public_name, donor_id, supplier_name, control, taxon_id, created from current_samples where name = ?];
            $sample_sth = $dbh->prepare($sql);
            $sample_sth->execute;
            $sample_sth->bind_columns(\($public_name, $donor_id, $supplier_name, $control, $taxon_id, $created));
            
            $sql       = q[select name from current_studies where internal_id = ?];
            $study_sth = $dbh->prepare($sql);
            $study_sth->execute;
            $study_sth->bind_col(1, \$study_title);
        }
        
        my $cmd = "imeta -z $zone qu -d $query";
        open(my $qu_fh, $cmd . ' |') || $self->throw("Could not open pipe to [$cmd]");
        my %files;
        my $collection;
        my (%analysis_to_col_date, %analysis_files);
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
                        if (exists $meta->{$attribute}) {
                            my $previous = $meta->{$attribute};
                            if (ref($previous)) {
                                $meta->{$attribute} = [@$previous, $1];
                            }
                            else {
                                $meta->{$attribute} = [$previous, $1];
                            }
                        }
                        else {
                            $meta->{$attribute} = $1;
                        }
                    }
                }
                close($ls_fh) || $self->throw("Could not close pipe to [$ls_cmd]");
                
                # grab extra info from the warehouse database if requested
                if ($sample_sth) {
                    undef $public_name;
                    undef $donor_id;
                    undef $supplier_name;
                    undef $control;
                    undef $taxon_id;
                    undef $created;
                    my $sanger_sample_id = $meta->{sample};
                    if ($sanger_sample_id) {
                        $sample_sth->execute($sanger_sample_id);
                        $sample_sth->fetch;
                        if ($public_name) {
                            $meta->{public_name} = "$public_name";
                        }
                        if ($donor_id) {
                            $meta->{sample_cohort} = "$donor_id";
                        }
                        if ($supplier_name) {
                            $meta->{sample_supplier_name} = "$supplier_name";
                        }
                        if (defined $control) {
                            $meta->{sample_control} = "$control";
                        }
                        if (defined $taxon_id) {
                            $meta->{taxon_id} = "$taxon_id";
                        }
                        if (defined $created) {
                            $meta->{sample_created_date} = "$created";
                        }
                    }
                    
                    my $study_id = $meta->{study_id};
                    if ($study_id) {
                        my @study_ids = ref($study_id) ? @$study_id : ($study_id);
                        
                        if (@study_ids > 1 && $query =~ /study_id = (\d+)/) {
                            # because we most likely want to populate VRTrack
                            # based on this metadata, and VRTrack can't cope
                            # with multiple study_ids per sample, check the
                            # source to see if we wanted just a single study and
                            # limit to that one
                            my $desired = $1;
                            foreach my $study_id (@study_ids) {
                                if ($study_id == $desired) {
                                    @study_ids = ($desired);
                                    $meta->{study_id} = $desired;
                                    last;
                                }
                            }
                        }
                        
                        my @study_titles;
                        foreach my $study_id (@study_ids) {
                            undef $study_title;
                            $study_sth->execute($study_id);
                            $study_sth->fetch;
                            if ($study_title) {
                                push(@study_titles, "$study_title");
                            }
                        }
                        
                        if (@study_titles > 1) {
                            $meta->{study_title} = \@study_titles;
                        }
                        else {
                            $meta->{study_title} = $study_titles[0];
                        }
                    }
                    
                    # another Sanger-specific thing is if we had an
                    # analysis_uuid, see if there are any collections with the
                    # same analysis_uuid and store some of the associated files
                    if (exists $meta->{analysis_uuid}) {
                        # get the most recent collection associated with an
                        # analysis done for this file
                        my %collections;
                        my @uuids = ref($meta->{analysis_uuid}) ? @{ $meta->{analysis_uuid} } : ($meta->{analysis_uuid});
                        foreach my $uuid (@uuids) {
                            my ($col, $created_date);
                            if (exists $analysis_to_col_date{$uuid}) {
                                ($col, $created_date) = @{ $analysis_to_col_date{$uuid} };
                            }
                            else {
                                my $colls_cmd = "imeta -z $zone qu -C analysis_uuid = $uuid";
                                open(my $colls_fh, $colls_cmd . ' |') || $self->throw("Could not open pipe to [$colls_cmd]");
                                while (<$colls_fh>) {
                                    if (/^collection:\s+(\S+)/) {
                                        my $this_col = $1;
                                        my $date_cmd = "imeta ls -C $this_col dcterms:created";
                                        open(my $date_fh, $date_cmd . ' |') || $self->throw("Could not open pipe to [$date_cmd]");
                                        while (<$date_fh>) {
                                            if (/^value: (.+)$/) {
                                                $created_date = $1;
                                            }
                                        }
                                        close($date_fh) || $self->throw("Could not close pipe to [$date_cmd]");
                                        $created_date   || next;
                                        $col = $this_col;
                                        #*** do we ever get multiple collections
                                        #    per uuid?
                                        last;
                                    }
                                }
                                close($colls_fh) || $self->throw("Could not close pipe to [$colls_cmd]");
                                
                                $col || next;
                                $analysis_to_col_date{$uuid} = [$col, $created_date];
                            }
                            
                            $collections{$col} = $created_date;
                        }
                        my ($analysis_collection) = sort { $collections{$b} cmp $collections{$a} } keys %collections;
                        
                        if ($analysis_collection) {
                            my @files;
                            if (exists $analysis_files{$analysis_collection}) {
                                @files = @{ $analysis_files{$analysis_collection} };
                            }
                            else {
                                my $ils_cmd = "ils -r $analysis_collection";
                                open(my $ils_fh, $ils_cmd . ' |') || $self->throw("Could not open pipe to [$ils_cmd]");
                                my $dir;
                                while (<$ils_fh>) {
                                    chomp;
                                    if (/^($analysis_collection[^:]*)/) {
                                        $dir = $1;
                                    }
                                    elsif (/(\S+(?:Sample_Probe_Profile|annotation|\.fcr\.)\S+)/i) {
                                        # we could store all files, or take the
                                        # above regex as a user-option, but for now
                                        # we'll just hard-code it since this is all
                                        # Sanger-specific anyway
                                        next if /quantile/;
                                        push(@files, file($dir, $1)->stringify);
                                    }
                                }
                                close($ils_fh) || $self->throw("Could not close pipe to [$ils_cmd]");
                                
                                $analysis_files{$analysis_collection} = [@files];
                            }
                            
                            if (@files) {
                                $meta->{irods_analysis_files} = @files > 1 ? \@files : $files[0];
                            }
                        }
                    }
                }
                
                $files{$path} = $meta;
            }
        }
        close($qu_fh) || $self->throw("Could not close pipe to [$cmd]");
        
        return \%files;
    }
    
    method all (Defined :$handle!, Str :$file_query!, Str|Dir :$local_root_dir!, Str :$update_interval?) {
        my %args;
        $args{handle}          = $handle;
        $args{file_query}      = $file_query;
        $args{local_root_dir}  = $local_root_dir;
        $args{update_interval} = $update_interval if defined($update_interval);
        return $self->_all_files(%args);
    }
    
    method all_with_warehouse_metadata (Defined :$handle!, Str :$file_query!, Str|Dir :$local_root_dir!, Str :$update_interval?) {
        my %args;
        $args{handle}          = $handle;
        $args{file_query}      = $file_query;
        $args{local_root_dir}  = $local_root_dir;
        $args{update_interval} = $update_interval if defined($update_interval);
        return $self->_all_files(%args, add_metadata_from_warehouse => 1);
    }
    
    method _all_files (Defined :$handle!, Str :$file_query!, Str|Dir :$local_root_dir!, Str :$update_interval?, Bool :$add_metadata_from_warehouse?) {
        # _get_irods_files_and_metadata will get called twice in row: once to
        # see if the datasource changed, and again here; _has_changed caches
        # the result, and we clear the cache after getting that data
        $add_metadata_from_warehouse ||= 0;
        my $files = $self->_get_irods_files_and_metadata($handle, $file_query, $add_metadata_from_warehouse);
        $self->_clear_cache;
        
        my %ignore_keys = map { $_ => 1 } qw(study_id study_title sample_common_name ebi_sub_acc reference ebi_sub_md5 ebi_run_acc ebi_sub_date manual_qc sample_created_date taxon_id);
        
        my $did = $self->_datasource_id;
        my @element_args;
        foreach my $path (sort { $files->{$a}->{vrpipe_irods_order} <=> $files->{$b}->{vrpipe_irods_order} } keys %$files) {
            my $new_metadata = $files->{$path};
            delete $new_metadata->{vrpipe_irods_order};
            
            my $sub_path = $path;
            $sub_path =~ s/^\///;
            my $file_abs_path = file($local_root_dir, $sub_path)->stringify;
            
            # consider type to be any if not defined in the irods metadata; if
            # not a VRPipe filetype it will be treated as an any
            my $type = delete $new_metadata->{type};
            $type ||= 'any';
            
            my $vrfile = VRPipe::File->create(path => $file_abs_path, type => $type)->original;
            
            # add metadata to file, detecting any changes
            my $current_metadata = $vrfile->metadata;
            my $changed          = 0;
            if ($current_metadata && keys %$current_metadata) {
                while (my ($key, $val) = each %$current_metadata) {
                    if ($add_metadata_from_warehouse) {
                        next if exists $ignore_keys{$key};
                    }
                    
                    next unless defined $val;
                    next unless defined $new_metadata->{$key};
                    if (_vals_different($val, $new_metadata->{$key})) {
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
    
    sub _vals_different {
        my ($orig, $new) = @_;
        if (!ref($orig) && !ref($new)) {
            return $orig ne $new;
        }
        
        if ((ref($orig) && !ref($new)) || (!ref($orig) && ref($new))) {
            return 1;
        }
        
        my %orig = map { $_ => 1 } @$orig;
        my %new  = map { $_ => 1 } @$new;
        foreach my $orig (keys %orig) {
            return 1 unless delete $new{$orig};
        }
        return 1 if keys %new;
        
        return 0;
    }
}

1;
