use VRPipe::Base;

class VRPipe::DataSource::vrtrack with VRPipe::DataSourceRole {
    use VertRes::Utils::VRTrackFactory;
    use VertRes::Utils::Hierarchy;
    
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
        #*** need to figure out a nice way of determining if anything has
        #    changed in the VRTrack db...
        return 1;
    }
    
    method _update_changed_marker {
        #*** this needs to be properly implemented as well
        $self->_changed_marker('changed');
    }
    
    method lanes (Defined :$handle,
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
        
        return \@elements;
    }
}

=pod
    =head2 next_lane_path
    
      Arg [1]    : hash of hierarchy level name keys, with regex values. Special
                   key of 'processed' with hashref value as per
                   processed_lane_hnames.
      Example    : all lanes:
                   while (my $hname = $track->next_lane_path()) { ... }
                   mapped, not yet improved lanes for sample NA01:
                   while (my $hname = $track->next_lane_hname(sample => 'NA01',
                       processed => { mapped => 1, improved => 0 })) { ... }
      Description: retrieves a (optionally filtered) list of all lane hierarchy
                   paths, ordered by project, sample, library names.
      Returntype : arrayref
    
    =cut
    
    sub next_lane_path {
        my ($self, %args) = @_;
        my $processed = delete $args{processed};
        my $store_name;
        foreach my $key (sort keys %args) {
            $store_name .= $key.'->'.$args{$key};
        }
        if ($processed) {
            foreach my $key (sort keys %{$processed}) {
                $store_name .= $key.'->'.$processed->{$key};
            }
        }
        
        unless (exists $self->{$store_name}) {
            
        }
        
        if (exists $self->{$store_name}) {
            my $sth = $self->{$store_name};
            my $row_data = $sth->fetchrow_arrayref;
            if ($row_data) {
                
            }
            else {
                delete $self->{$store_name};
                return;
            }
        }
    }
=cut

1;