
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
    use VRPipe::Schema;
    use VRPipe::Steps::irods;
    
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
            return "An element will comprise one of the files returned by imeta qu -d given the arguments you supply for the 'file_query' option (which can be the options specified directly, or the absolute path to a file containing multiple sets of imeta options, 1 set per line). The file will have all the relevant irods metadata associated with it, and a local path based on the 'local_root_dir' option. To avoid spamming the irods server, the update_interval option allows you to specify the minimum number of minutes between each check for changes to files. If update_interval is not supplied it defaults to 5 seconds when testing, and 1 day in production.";
        }
        if ($method eq 'all_with_warehouse_metadata') {
            return "In addition to doing everything the all method does, it adds extra metadata found in the warehouse database to the files with the keys public_name, sample_supplier_name, sample_control, sample_cohort, taxon_id, sample_created_date and study_title (if defined). Optionally provide a comma-separated list of required keys to the required_metadata option to ignore files lacking that metadata. If any analysis has been done to a file, the associated files are stored under irods_analysis_files. The vrtrack_group option determines which group the studies your data are in are placed under - use it for grouping together studies you will analyse the same way later. (This method is Sanger-specific and also requires the environment variables WAREHOUSE_DATABASE, WAREHOUSE_HOST, WAREHOUSE_PORT and WAREHOUSE_USER.)";
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
        my $update_interval = $options->{update_interval} || $VRPipe::Persistent::InMemory::deployment eq 'production' ? 1440 : 999999; # in mins
        $update_interval *= 60;                                                                                                         # in seconds
        my $lock_key         = 'irods_datasource.' . $self->_datasource_id;
        my $locked           = $im->lock($lock_key, unlock_after => $update_interval);
        my $current_checksum = $self->_changed_marker;
        
        if ($current_checksum && length($current_checksum) == 32) {
            return $current_checksum unless $locked;
        }
        # else we always get the latest checksum if we have no valid checksum
        
        # get the current files and their metadata and stringify it all
        my $files = $self->_get_irods_files_and_metadata($self->_open_source(), $options->{file_query}, $self->method eq 'all_with_warehouse_metadata', $options->{required_metadata}, $options->{vrtrack_group});
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
    
    method _get_irods_files_and_metadata (Str $zone!, Str $raw_query!, Bool $add_metadata_from_warehouse?, Maybe[Str] $required_metadata?, Maybe[Str] $vrtrack_group?) {
        return $self->_irods_files_and_metadata_cache if $self->_cached;
        my @required_keys;
        if ($required_metadata) {
            @required_keys = split(',', $required_metadata);
        }
        
        my ($sample_sth, $public_name, $donor_id, $supplier_name, $control, $taxon_id, $created, $gender);
        my ($study_sth, $study_title);
        my $vrtrack;
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
            
            my $sql = q[select public_name, donor_id, supplier_name, control, taxon_id, created, gender from current_samples where name = ?];
            $sample_sth = $dbh->prepare($sql);
            $sample_sth->execute;
            $sample_sth->bind_columns(\($public_name, $donor_id, $supplier_name, $control, $taxon_id, $created, $gender));
            
            $sql       = q[select name from current_studies where internal_id = ?];
            $study_sth = $dbh->prepare($sql);
            $study_sth->execute;
            $study_sth->bind_col(1, \$study_title);
            
            $vrtrack = VRPipe::Schema->create('VRTrack');
            $vrtrack_group ||= 'all_studies';
        }
        
        my %files;
        my @queries;
        if ($raw_query !~ /\s/ && -e $raw_query) {
            # it's a file containing multiple queries
            my $file = VRPipe::File->create(path => $raw_query);
            my $fh = $file->openr;
            while (<$fh>) {
                chomp;
                push(@queries, $_);
            }
            $file->close;
        }
        else {
            # it's a single directly-specified query
            @queries = ($raw_query);
        }
        
        my (%analysis_to_cols, %col_dates, %analysis_files, %analyses, %collections);
        my $order = 1;
        foreach my $query (@queries) {
            my @cmd_output = VRPipe::Steps::irods->open_irods_command("imeta -z $zone qu -d $query");
            my $collection;
            QU: foreach (@cmd_output) {
                #*** do we have to worry about spaces in file paths?...
                if (/^collection:\s+(\S+)/) {
                    $collection = $1;
                }
                elsif (/^----/) {
                    undef $collection;
                }
                elsif ($collection && /^dataObj:\s+(\S+)/) {
                    my $path = "$collection/$1";
                    # get all the metadata for this file
                    my $meta = VRPipe::Steps::irods->get_file_metadata($path);
                    $meta->{vrpipe_irods_order} = $order++;
                    
                    # grab extra info from the warehouse database if requested
                    if ($sample_sth) {
                        undef $public_name;
                        undef $donor_id;
                        undef $supplier_name;
                        undef $control;
                        undef $taxon_id;
                        undef $created;
                        undef $gender;
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
                            if (defined $gender) {
                                $meta->{sample_gender} = "$gender";
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
                            my @uuids = ref($meta->{analysis_uuid}) ? @{ $meta->{analysis_uuid} } : ($meta->{analysis_uuid});
                            my %collections;
                            foreach my $uuid (@uuids) {
                                if (exists $analysis_to_cols{$uuid}) {
                                    foreach my $this_col (@{ $analysis_to_cols{$uuid} }) {
                                        $collections{$this_col} = 1;
                                    }
                                }
                                else {
                                    my $vrtrack_analysis = $vrtrack->add('Analysis', { uuid => $uuid });
                                    $analyses{$uuid} = $vrtrack_analysis;
                                    my @cmd_output = VRPipe::Steps::irods->open_irods_command("imeta -z $zone qu -C analysis_uuid = $uuid");
                                    
                                    my @these_cols;
                                    foreach (@cmd_output) {
                                        if (/^collection:\s+(\S+)/) {
                                            my $this_col = $1;
                                            $collections{$this_col} = 1;
                                            push(@these_cols, $this_col);
                                            
                                            unless (exists $col_dates{$this_col}) {
                                                my @date_cmd_output = VRPipe::Steps::irods->open_irods_command("imeta ls -C $this_col dcterms:created");
                                                my $date;
                                                foreach (@date_cmd_output) {
                                                    if (/^value: (.+)$/) {
                                                        $date = $1;
                                                        last;
                                                    }
                                                }
                                                $date ||= '2013-01-01T12:00:00';
                                                $col_dates{$this_col} = $date;
                                                
                                                my $vrtrack_col = $vrtrack->add('Collection', { path => $this_col, date => $vrtrack->date_to_epoch($date) }, incoming => { type => 'directory', node => $vrtrack_analysis });
                                                $collections{$this_col} = $vrtrack_col;
                                            }
                                        }
                                    }
                                    
                                    $analysis_to_cols{$uuid} = \@these_cols;
                                }
                            }
                            
                            my ($analysis_collection) = sort { $col_dates{$b} cmp $col_dates{$a} } keys %collections;
                            
                            if ($analysis_collection) {
                                my @files;
                                if (exists $analysis_files{$analysis_collection}) {
                                    @files = @{ $analysis_files{$analysis_collection} };
                                }
                                else {
                                    my @cmd_output = VRPipe::Steps::irods->open_irods_command("ils -r $analysis_collection");
                                    my $dir;
                                    foreach (@cmd_output) {
                                        if (/^($analysis_collection[^:]*)/) {
                                            $dir = $1;
                                        }
                                        elsif (/^\s+(\w\S+)$/) {
                                            #*** there are files with spaces in
                                            # the filename and also ones that
                                            # start with special chars like ~$,
                                            # but it's too awkward to bother
                                            # supporting them - they are ignored!
                                            my $path = file($dir, $1)->stringify;
                                            push(@files, $path);
                                            
                                            my $vrtrack_file = $vrtrack->add_file($path);
                                            $collections{$analysis_collection}->relate_to($vrtrack_file, 'contains');
                                        }
                                    }
                                    
                                    $analysis_files{$analysis_collection} = [@files];
                                }
                                
                                if (@files) {
                                    $meta->{irods_analysis_files} = @files > 1 ? \@files : $files[0];
                                }
                            }
                        }
                    }
                    
                    foreach my $key (@required_keys) {
                        next QU unless defined $meta->{$key};
                    }
                    
                    # represent lab tracking metadata in the graph database
                    if ($vrtrack) {
                        # general
                        my $group = $vrtrack->add('Group', { name => $vrtrack_group });
                        my $study = $vrtrack->add('Study', { id => $meta->{study_id}, name => $meta->{study_title} || $meta->{study}, accession => $meta->{study_accession_number} }, incoming => { type => 'has', node => $group });
                        
                        my $donor;
                        if (defined $meta->{sample_cohort}) {
                            $donor = $vrtrack->add('Donor', { id => $meta->{sample_cohort} }, incoming => { type => 'member', node => $study });
                        }
                        
                        my $sample_created_date;
                        if (defined $meta->{sample_created_date}) {
                            # convert '2013-05-10 06:45:32' to epoch seconds
                            $sample_created_date = $vrtrack->date_to_epoch($meta->{sample_created_date});
                        }
                        my $sample = $vrtrack->add('Sample', { name => $meta->{sample}, public_name => $meta->{public_name}, id => $meta->{sample_id}, supplier_name => $meta->{sample_supplier_name}, accession => $meta->{sample_accession_number}, created_date => $sample_created_date, consent => $meta->{sample_consent}, control => $meta->{sample_control} }, incoming => { type => 'member', node => $study });
                        
                        my $sex = 'U'; # unknown
                        if (defined $meta->{sample_gender}) {
                            if ($meta->{sample_gender} =~ /^f/i) {
                                $sex = 'F';
                            }
                            elsif ($meta->{sample_gender} =~ /^m/i) {
                                $sex = 'M';
                            }
                            delete $meta->{sample_gender}; # this new thing we don't want on our VRPipe::File metadata
                        }
                        $vrtrack->add('Gender', { source_gender_md5 => md5_hex('sequencescape.' . $sex), source => 'sequencescape', gender => $sex }, incoming => { type => 'gender', node => $sample });
                        
                        if ($donor) {
                            #*** we don't have a way of tracking changes to
                            # relationships, so won't have a record if we swap
                            # to a different donor here
                            $donor->relate_to($sample, 'sample', selfish => 1);
                        }
                        
                        if (defined $meta->{taxon_id}) {
                            my $taxon = $vrtrack->add('Taxon', { id => $meta->{taxon_id}, common_name => $meta->{sample_common_name} });
                            $taxon->relate_to($sample, 'member', selfish => 1);
                        }
                        
                        my $file = $vrtrack->add_file($path);
                        $file->add_properties({ manual_qc => $meta->{manual_qc}, target => $meta->{target}, md5 => $meta->{md5} });
                        my $file_connected = 0;
                        my $unique         = file($path)->basename;
                        $unique =~ s/\.gz$//;
                        $unique =~ s/\.[^\.]+$//;
                        
                        if (defined $meta->{analysis_uuid}) {
                            foreach my $uuid (ref($meta->{analysis_uuid}) ? @{ $meta->{analysis_uuid} } : ($meta->{analysis_uuid})) {
                                $file->relate_to($analyses{$uuid}, 'analysed');
                            }
                        }
                        
                        # bams
                        my $library;
                        if (defined $meta->{library_id}) {
                            $library = $vrtrack->add('Library', { id => $meta->{library_id}, name => $meta->{library}, tag => $meta->{tag} });
                            $sample->relate_to($library, 'prepared', selfish => 1);
                        }
                        
                        my $lane;
                        if (defined $meta->{lane}) {
                            $lane = $vrtrack->add('Lane', { unique => $unique, lane => $meta->{lane}, run => $meta->{id_run}, total_reads => $meta->{total_reads}, is_paired_read => $meta->{is_paired_read} });
                            
                            $lane->relate_to($file, 'aligned', selfish => 1);
                            $file_connected = 1;
                            
                            if ($library) {
                                $library->relate_to($lane, 'sequenced', selfish => 1);
                            }
                        }
                        
                        if ($meta->{alignment} && defined $meta->{reference}) {
                            my $aln = $vrtrack->add('Alignment', { reference => $meta->{reference} });
                            $aln->relate_to($file, 'reference', selfish => 1);
                        }
                        
                        if (defined $meta->{ebi_sub_acc} && $meta->{ebi_run_acc}) {
                            my $sub = $vrtrack->add('EBI_Submission', { acc => $meta->{ebi_sub_acc}, defined $meta->{ebi_sub_date} ? (sub_date => $vrtrack->date_to_epoch($meta->{ebi_sub_date})) : () });
                            
                            my $run = $vrtrack->add('EBI_Run', { acc => $meta->{ebi_run_acc}, md5 => $meta->{ebi_sub_md5} }, incoming => { type => 'submitted', node => $file }, outgoing => { type => 'submitted', node => $sub });
                        }
                        
                        # infinium idats and gtc
                        if (defined $meta->{beadchip}) {
                            my $beadchip = $vrtrack->add('Beadchip', { id => $meta->{beadchip}, design => $meta->{beadchip_design} }, incoming => { type => 'has', node => $study });
                            
                            if (defined $meta->{beadchip_section}) {
                                my $section = $vrtrack->add('Section', { unique => $unique, section => $meta->{beadchip_section} }, incoming => { type => 'placed', node => $sample });
                                $beadchip->relate_to($section, 'section', selfish => 1);
                                $section->relate_to($file, 'processed', selfish => 1);
                                $file_connected = 1;
                            }
                        }
                        
                        my $isample;
                        if (defined $meta->{infinium_sample}) {
                            $isample = $vrtrack->add('Infinium_Sample', { id => $meta->{infinium_sample} });
                            $sample->relate_to($isample, 'infinium', selfish => 1);
                        }
                        
                        if (defined $meta->{infinium_plate}) {
                            my $plate = $vrtrack->add('Infinium_Plate', { id => $meta->{infinium_plate} }, incoming => { type => 'has', node => $study });
                            
                            if (defined $meta->{infinium_well}) {
                                my $well = $vrtrack->add('Well', { unique => $meta->{infinium_plate} . '.' . $meta->{infinium_well}, well => $meta->{infinium_well} }, incoming => { type => 'placed', node => $isample || $sample });
                                $plate->relate_to($well, 'well', selfish => 1);
                            }
                        }
                        
                        unless ($file_connected) {
                            $sample->relate_to($file, 'processed');
                        }
                    }
                    
                    $files{$path} = $meta;
                }
            }
        }
        
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
    
    method all_with_warehouse_metadata (Defined :$handle!, Str :$file_query!, Str|Dir :$local_root_dir!, Str :$update_interval?, Str :$required_metadata?, Str :$vrtrack_group?) {
        my %args;
        $args{handle}          = $handle;
        $args{file_query}      = $file_query;
        $args{local_root_dir}  = $local_root_dir;
        $args{update_interval} = $update_interval if defined($update_interval);
        return $self->_all_files(%args, add_metadata_from_warehouse => 1, $required_metadata ? (required_metadata => $required_metadata) : (), $vrtrack_group ? (vrtrack_group => $vrtrack_group) : ());
    }
    
    method _all_files (Defined :$handle!, Str :$file_query!, Str|Dir :$local_root_dir!, Str :$update_interval?, Bool :$add_metadata_from_warehouse?, Str :$required_metadata?, Str :$vrtrack_group?) {
        # _get_irods_files_and_metadata will get called twice in row: once to
        # see if the datasource changed, and again here; _has_changed caches
        # the result, and we clear the cache after getting that data
        $add_metadata_from_warehouse ||= 0;
        my $files = $self->_get_irods_files_and_metadata($handle, $file_query, $add_metadata_from_warehouse, $required_metadata, $vrtrack_group);
        $self->_clear_cache;
        
        my %ignore_keys = map { $_ => 1 } qw(study_id study_title sample_common_name ebi_sub_acc reference ebi_sub_md5 ebi_run_acc ebi_sub_date sample_created_date taxon_id lane);
        
        my $did = $self->_datasource_id;
        my @element_args;
        my @changed_details;
        my $anti_repeat_store = {};
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
                    if (my $diff = $self->_vals_different($val, $new_metadata->{$key})) {
                        $changed = 1;
                        push(@changed_details, "$file_abs_path $key: $diff");
                        last;
                    }
                }
            }
            if ($add_metadata_from_warehouse && exists $new_metadata->{lane}) {
                my $lane = $vrfile->basename;
                $lane =~ s/\.gz$//;
                $lane =~ s/\.[^\.]+$//;
                $new_metadata->{lane} = $lane;
            }
            
            # if there was no metadata this will add metadata to the file.
            $vrfile->add_metadata($new_metadata, replace_data => 0);
            $vrfile->add_metadata({ irods_path => $path, defined $new_metadata->{md5} ? (expected_md5 => $new_metadata->{md5}) : () });
            
            my $result_hash = { paths => [$file_abs_path], irods_path => $path };
            if ($changed) {
                $result_hash->{changed} = [[$vrfile, $new_metadata]];
                $self->_start_over_elements_due_to_file_metadata_change($result_hash, \@changed_details, $anti_repeat_store);
                delete $result_hash->{changed};
            }
            push(@element_args, { datasource => $did, result => $result_hash });
        }
        $self->_create_elements(\@element_args);
    }
}

1;
