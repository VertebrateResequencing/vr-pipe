use VRPipe::Base;

class VRPipe::DataSource::vrtrack with VRPipe::DataSourceRole {
    # eval these so that test suite can pass syntax check on this module when
    # VertRes is not installed
    eval "use VertRes::Utils::VRTrackFactory;";
    eval "use VertRes::Utils::Hierarchy;";
    use Digest::MD5 qw(md5_hex);
    
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
      my $lane_change = VRTrack::Lane->_all_values_by_field($vrtrack_source, 'changed');
      my $file_md5    = VRTrack::File->_all_values_by_field($vrtrack_source, 'md5');
      my $digest      = md5_hex join( @$lane_change, map { defined $_ ? $_ : 'NULL' } @$file_md5); 
      return $digest;
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
        
        my @elements;
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
            
            push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => {lane => $lane->hierarchy_name}, withdrawn => 0));
        }
        $self->_update_changed_marker; 
        return \@elements;
    }

    method lane_fastqs (Defined :$handle!,
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
                  Bool :$import?,
                  Bool :$qc?,
                  Bool :$mapped?,
                  Bool :$stored?,
                  Bool :$deleted?,
                  Bool :$swapped?,
                  Bool :$altered_fastq?,
                  Bool :$improved?,
                  Bool :$snp_called?) {
        #push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => {paths => $hash_ref->{paths}, lane => $lane}, withdrawn => 0));
   }
}
1;
