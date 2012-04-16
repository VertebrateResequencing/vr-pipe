use VRPipe::Base;

class VRPipe::DataSource::vrtrack with VRPipe::DataSourceRole {
    # eval these so that test suite can pass syntax check on this module when
    # VertRes is not installed
    eval "use VertRes::Utils::VRTrackFactory;";
    eval "use VertRes::Utils::Hierarchy;";
    use Digest::MD5 qw(md5_hex);
    use File::Spec::Functions;
    
    our %file_type_to_type = (0 => 'fq', 1 => 'fq', 2 => 'fq', 3 => 'fq', 4 => 'bam', 5 => 'bam', 6 => 'cram');
    
    method description {
        return "Use a VRTrack database to extract information from";
    }
    method source_description {
        return "The name of the VRTrack database; assumes your database connection details are held in the normal set of VRTrack-related environment variables";
    }
    method method_description (Str $method) {
        if ($method eq 'lanes') {
            return "An element will comprise the name of a lane (only).";
        }
        elsif ($method eq 'lane_bams') {
            return "An element will comprise all the bams for a single lane, and the bam files will have all relevant available metadata associated with them.";
        }
        elsif ($method eq 'lane_fastqs') {
            return "An element will comprise all the fastqs for a single lane, and the fastq files will have all relevant available metadata associated with them.";
        }
        
        return '';
    }
    
    method _open_source {
        return VertRes::Utils::VRTrackFactory->instantiate(database => $self->source, mode => 'r');
    }
    
    method _has_changed {
        return 1 unless defined($self->_changed_marker);#on first instantiation _changed_marker is undefined, defaults to changed in this case 
        return 1 if ($self->_vrtrack_lane_file_checksum ne $self->_changed_marker);#checks for new or deleted lanes or changed files(including deleted/added files)
        return 0; 
    }
    
    method _update_changed_marker { 
        $self->_changed_marker($self->_vrtrack_lane_file_checksum); 
    }

    method _vrtrack_lane_file_checksum {
        my $vrtrack_source = $self->_open_source();
        
        # concatenating the 'changed' column of all lanes would let us know if
        # anything at all changed in the lanes, but means that if we're running
        # a pipeline that changes something in the lane table we'll have to
        # rebuild our dataelemnts constantly, which is prohibitvely too slow.
        # Instead we look at the lane-related options the user will be filtering
        # by and only check to see if those related columns have changed. To
        # detect overall changes to lanes (gain/loss) we always also check
        # lane_id and withdrawn. Things like raw_reads changing are irrelvant;
        # important changes to the actual file data will be caught by the file
        # md5s changing
        
        # we always care about lane ids changing
        my $lane_changes = VRTrack::Lane->_all_values_by_field($vrtrack_source, 'lane_id', 'hierarchy_name');
        push(@$lane_changes, @{VRTrack::Lane->_all_values_by_field($vrtrack_source, 'withdrawn', 'hierarchy_name')});
        
        # we care about options if supplied
        my $options = $self->options;
        foreach my $opt (qw(import qc mapped stored deleted swapped altered_fastq improved snp_called)) {
            if (defined $options->{$opt}) {
                push(@$lane_changes, @{VRTrack::Lane->_all_values_by_field($vrtrack_source, 'processed', 'hierarchy_name')});
                last;
            }
        }
        push(@$lane_changes, @{VRTrack::Lane->_all_values_by_field($vrtrack_source, 'qc_status', 'hierarchy_name')}) if defined $options->{qc_status};
        push(@$lane_changes, @{VRTrack::Lane->_all_values_by_field($vrtrack_source, 'gt_status', 'hierarchy_name')}) if defined $options->{gt_status};
        
        # we care about certain types of files depending on method
        my $method = $self->method;
        if ($method eq 'lane_fastqs') {
            push(@$lane_changes, @{VRTrack::File->_all_values_by_field($vrtrack_source, 'md5', 'hierarchy_name', '(type=0 or type=1 or type=2) and latest=true')});
        }
        elsif ($method eq 'lane_bams') {
            push(@$lane_changes, @{VRTrack::File->_all_values_by_field($vrtrack_source, 'md5', 'hierarchy_name', 'type=4 and latest=true')});
        }
        
        my $digest = md5_hex join('', map { defined $_ ? $_ : 'NULL' } @$lane_changes);
        return $digest;
    }
    
    method _filtered_lanes (Defined :$handle!,
                            ArrayRef :$project?,
                            ArrayRef :$sample?,
                            ArrayRef :$individual?,
                            ArrayRef :$population?,
                            ArrayRef :$platform?,
                            ArrayRef :$centre?,
                            ArrayRef :$library?,
                            Str :$project_regex?,
                            Str :$sample_regex?,
                            Str :$library_regex?,
                            Str :$qc_status?,
                            Str :$gt_status?,
                            Bool :$import?,
                            Bool :$qc?,
                            Bool :$mapped?,
                            Bool :$stored?,
                            Bool :$deleted?,
                            Bool :$swapped?,
                            Bool :$altered_fastq?,
                            Bool :$improved?,
                            Bool :$snp_called?) {
        my $hu = VertRes::Utils::Hierarchy->new();
        my @lanes = $hu->get_lanes(vrtrack => $handle,
                                   $project ? (project => $project) : (),
                                   $sample ? (sample => $sample) : (),
                                   $individual ? (individual => $individual) : (),
                                   $population ? (population => $population) : (),
                                   $platform ? (platform => $platform) : (),
                                   $centre ? (centre => $centre) : (),
                                   $library ? (library => $library) : (),
                                   $project_regex ? (project_regex => $project_regex) : (),
                                   $sample_regex ? (sample_regex => $sample_regex) : (),
                                   $library_regex ? (library_regex => $library_regex) : ());
        
        my @filtered;
        foreach my $lane (@lanes) {
            if (defined $import) {
                my $processed = $lane->is_processed('import');
                next if $processed != $import;
            }
            if (defined $qc) {
                my $processed = $lane->is_processed('qc');
                next if $processed != $qc;
            }
            if (defined $mapped) {
                my $processed = $lane->is_processed('mapped');
                next if $processed != $mapped;
            }
            if (defined $stored) {
                my $processed = $lane->is_processed('stored');
                next if $processed != $stored;
            }
            if (defined $deleted) {
                my $processed = $lane->is_processed('deleted');
                next if $processed != $deleted;
            }
            if (defined $swapped) {
                my $processed = $lane->is_processed('swapped');
                next if $processed != $swapped;
            }
            if (defined $altered_fastq) {
                my $processed = $lane->is_processed('altered_fastq');
                next if $processed != $altered_fastq;
            }
            if (defined $improved) {
                my $processed = $lane->is_processed('improved');
                next if $processed != $improved;
            }
            if (defined $snp_called) {
                my $processed = $lane->is_processed('snp_called');
                next if $processed != $snp_called;
            }
            if (defined $qc_status) {
                my $this_status = $lane->qc_status;
                next if $this_status !~ /$qc_status/;
            }
            if (defined $gt_status) {
                my $this_status = $lane->gt_status;
                next if $this_status !~ /$gt_status/;
            }
            
            push(@filtered, $lane);
        }
        
        return @filtered;
    }
    
    method lanes (Defined :$handle!,
                  ArrayRef :$project?,
                  ArrayRef :$sample?,
                  ArrayRef :$individual?,
                  ArrayRef :$population?,
                  ArrayRef :$platform?,
                  ArrayRef :$centre?,
                  ArrayRef :$library?,
                  Str :$project_regex?,
                  Str :$sample_regex?,
                  Str :$library_regex?,
                  Str :$qc_status?,
                  Str :$gt_status?,
                  Bool :$import?,
                  Bool :$qc?,
                  Bool :$mapped?,
                  Bool :$stored?,
                  Bool :$deleted?,
                  Bool :$swapped?,
                  Bool :$altered_fastq?,
                  Bool :$improved?,
                  Bool :$snp_called?) {
        my %args;
        $args{handle} = $handle if defined($handle);
        $args{project} = $project if defined($project);
        $args{sample} = $sample if defined($sample);
        $args{individual} = $population if defined($population);
        $args{platform} = $platform if defined($platform);
        $args{centre} = $centre if defined($centre);
        $args{library} = $library if defined($library);
        $args{project_regex} = $project_regex if defined($project_regex);
        $args{sample_regex} = $sample_regex if defined($sample_regex);
        $args{library_regex} = $library_regex if defined($library_regex);
        $args{import} = $import if defined($import);
        $args{qc} = $qc if defined($qc);
        $args{mapped} = $mapped if defined($mapped);
        $args{stored} = $stored if defined($stored);
        $args{deleted} = $deleted if defined($deleted);
        $args{swapped} = $swapped if defined($swapped);
        $args{altered_fastq} = $altered_fastq if defined($altered_fastq);
        $args{improved} = $improved if defined($improved);
        $args{snp_called} = $snp_called if defined($snp_called);
        $args{qc_status} = $qc_status if defined($qc_status);
        $args{gt_status} = $gt_status if defined($gt_status);
        
        my @elements;
        foreach my $lane ($self->_filtered_lanes(%args)) {
            push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => {lane => $lane->hierarchy_name}, withdrawn => 0));
        }
        $self->_update_changed_marker; 
        return \@elements;
    }

    method lane_bams (Defined :$handle!,
                      Str|Dir :$local_root_dir!,
                      ArrayRef :$project?,
                      ArrayRef :$sample?,
                      ArrayRef :$individual?,
                      ArrayRef :$population?,
                      ArrayRef :$platform?,
                      ArrayRef :$centre?,
                      ArrayRef :$library?,
                      Str :$project_regex?,
                      Str :$sample_regex?,
                      Str :$library_regex?,
                      Str :$qc_status?,
                      Str :$gt_status?,
                      Bool :$import?,
                      Bool :$qc?,
                      Bool :$mapped?,
                      Bool :$stored?,
                      Bool :$deleted?,
                      Bool :$swapped?,
                      Bool :$altered_fastq?,
                      Bool :$improved?,
                      Bool :$snp_called?) {
        my %args;
        $args{handle} = $handle if defined($handle);
        $args{local_root_dir} = $local_root_dir if defined($local_root_dir);
        $args{project} = $project if defined($project);
        $args{sample} = $sample if defined($sample);
        $args{individual} = $population if defined($population);
        $args{platform} = $platform if defined($platform);
        $args{centre} = $centre if defined($centre);
        $args{library} = $library if defined($library);
        $args{project_regex} = $project_regex if defined($project_regex);
        $args{sample_regex} = $sample_regex if defined($sample_regex);
        $args{library_regex} = $library_regex if defined($library_regex);
        $args{import} = $import if defined($import);
        $args{qc} = $qc if defined($qc);
        $args{mapped} = $mapped if defined($mapped);
        $args{stored} = $stored if defined($stored);
        $args{deleted} = $deleted if defined($deleted);
        $args{swapped} = $swapped if defined($swapped);
        $args{altered_fastq} = $altered_fastq if defined($altered_fastq);
        $args{improved} = $improved if defined($improved);
        $args{snp_called} = $snp_called if defined($snp_called);
        $args{qc_status} = $qc_status if defined($qc_status);
        $args{gt_status} = $gt_status if defined($gt_status);
        
        # add to the argument list to filter on bam files
        $args{'file_type'} = 4;
        return $self->_lane_files(%args);
    }
    
    method lane_fastqs(Defined :$handle!,
                       Str|Dir :$local_root_dir!,
                       ArrayRef :$project?,
                       ArrayRef :$sample?,
                       ArrayRef :$individual?,
                       ArrayRef :$population?,
                       ArrayRef :$platform?,
                       ArrayRef :$centre?,
                       ArrayRef :$library?,
                       Str :$project_regex?,
                       Str :$sample_regex?,
                       Str :$library_regex?,
                       Str :$qc_status?,
                       Str :$gt_status?,
                       Bool :$import?,
                       Bool :$qc?,
                       Bool :$mapped?,
                       Bool :$stored?,
                       Bool :$deleted?,
                       Bool :$swapped?,
                       Bool :$altered_fastq?,
                       Bool :$improved?,
                       Bool :$snp_called?) {
        my %args;
        $args{handle} = $handle if defined($handle);
        $args{local_root_dir} = $local_root_dir if defined($local_root_dir);
        $args{project} = $project if defined($project);
        $args{sample} = $sample if defined($sample);
        $args{individual} = $population if defined($population);
        $args{platform} = $platform if defined($platform);
        $args{centre} = $centre if defined($centre);
        $args{library} = $library if defined($library);
        $args{project_regex} = $project_regex if defined($project_regex);
        $args{sample_regex} = $sample_regex if defined($sample_regex);
        $args{library_regex} = $library_regex if defined($library_regex);
        $args{import} = $import if defined($import);
        $args{qc} = $qc if defined($qc);
        $args{mapped} = $mapped if defined($mapped);
        $args{stored} = $stored if defined($stored);
        $args{deleted} = $deleted if defined($deleted);
        $args{swapped} = $swapped if defined($swapped);
        $args{altered_fastq} = $altered_fastq if defined($altered_fastq);
        $args{improved} = $improved if defined($improved);
        $args{snp_called} = $snp_called if defined($snp_called);
        $args{qc_status} = $qc_status if defined($qc_status);
        $args{gt_status} = $gt_status if defined($gt_status);
        
        # add to the argument list to filter on fastq files
        $args{'file_type'} = '0|1|2';
        return $self->_lane_files(%args);
    }
   
    method _lane_files {
        my (undef, %args) = @_;
        my $file_type = delete $args{file_type};
        unless (defined $file_type) {
            $self->throw("file_type is required");
        }
        my ($file_type_key) = split('|', $file_type);
        my $vrpipe_filetype = $file_type_to_type{$file_type_key};
        
        my $local_root_dir = delete $args{local_root_dir};
        unless (defined $local_root_dir) {
            $self->throw("local_root_dir is required");
        }
        
        my $vrtrack = $args{handle};
        
        my @elements;
        my $lane_changed_hash;
        my $hu = VertRes::Utils::Hierarchy->new();
        foreach my $lane ($self->_filtered_lanes(%args)) {
            my %lane_info = $hu->lane_info($lane->name, vrtrack => $vrtrack, no_coverage => 1);
            
            my @files;
            foreach my $file (@{$lane->files}) {
                next unless $file->type =~ /^($file_type)$/;
                
                my $file_abs_path = file($local_root_dir, $file->name)->stringify; 
                my $new_metadata = {
                    expected_md5 => $file->md5,
                    center_name => $lane_info{'centre'},
                    project => $lane_info{'project'},
                    study => $lane_info{'study'},
                    $lane_info{'species'} ? (species => $lane_info{'species'}) : (),
                    population => $lane_info{ 'population'},
                    individual => $lane_info{'individual'},
                    sample => $lane_info{'sample'},
                    platform => $lane_info{'seq_tech'},
                    library => $lane_info{ 'library' },
                    lane => $lane_info{'lane'},
                    withdrawn => $lane_info{ 'withdrawn' } || 0, #*** we don't actually handle withdrawn files properly atm; if all withdrawn we shouldn't create the element...
                    insert_size => $lane_info{'insert_size'} || 0,
                    reads => $file->raw_reads || 0, 
                    bases => $file->raw_bases || 0, 
                    paired => $lane_info{'vrlane'}->is_paired,
                    lane_id => $file->lane_id
                };
                
                my $vrfile = VRPipe::File->get(path => $file_abs_path, type => $vrpipe_filetype);
                
                # add metadata to file but ensure that we update any fields in
                # the new metadata
                my $current_metadata = $vrfile->metadata;
                my $changed = 0;
                if ($current_metadata && keys %$current_metadata) {
                    foreach my $meta (qw(expected_md5 lane project study center_name sample population platform library insert_size)) {
                        next unless $new_metadata->{$meta};
                        if (defined $current_metadata->{$meta} && $current_metadata->{$meta} ne $new_metadata->{$meta}) {
                            $self->debug("metadata '$meta' changed from $current_metadata->{$meta} to $new_metadata->{$meta} for file $file_abs_path, so will mark lane ".$lane->name." as changed");
                            $changed = 1;
                            last;
                        }
                    }
                }
                
                # if there was no metadata this will add metadata to the file.
                $vrfile->add_metadata($new_metadata, replace_data => 0); 
                
                # if there was a change in VRPipe::File metadata store it in a hash and 
                # change the metadata in the VRPipe::File later when more appropriate, 
                # having made sure the DataElement element states have been reset (see below)                            
                if ($changed) {
                    push (@{$lane_changed_hash->{$lane->id}}, [$vrfile,$new_metadata]);
                }   
                
                push @files, $file_abs_path;             
            }
            
            push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result=>{ paths => \@files, lane=> $lane->name }, withdrawn => 0));
            
            if ($lane_changed_hash->{$lane->id}) {
                # reset element states first
                foreach my $estate ($elements[-1]->element_states) {
                    $estate->start_from_scratch;
                }
                # then change metadata in files
                foreach my $fm (@{$lane_changed_hash->{$lane->id}}) {
                    my ($vrfile, $new_metadata) = @$fm;
                    $vrfile->add_metadata($new_metadata, replace_data => 1);
                }
            }    
        } 
        return \@elements;
    }
}

1;