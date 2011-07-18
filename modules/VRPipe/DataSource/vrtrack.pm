use VRPipe::Base;

class VRPipe::DataSource::vrtrack with VRPipe::DataSourceRole {
    use VertRes::Utils::VRTrackFactory;
    use VertRes::Utils::Hierarchy;
    
    has '_hu_store' => (is => 'rw',
                        isa => 'ArrayRef',
                        predicate => '_hu_stored',
                        clearer => '_reset_hu',
                        traits  => ['Array'],
                        handles => { '_hu_push' => 'push',
                                     '_hu_shift' => 'shift' });
    
    method _open_source {
        return VertRes::Utils::VRTrackFactory->instantiate(database => $self->source, mode => 'r');
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
                  Bool :$snp_called?,
                  Bool :$return_objects = 0) {
        my $result;
        
        #*** need a proper, fast, efficient next_lane_path() method in VRTrack
        #    that will filter on hierarchy level and processed flag, and will
        #    access a single row from the db with each call...
        #    In the mean time just call HU get_lanes, manually filter results,
        #    and store array ref of results
        if ($self->_hu_stored) {
            $result = $self->_hu_shift;
        }
        else {
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
                
                $self->_hu_push($return_objects ? $lane : $lane->hierarchy_name);
            }
            
            $result = $self->_hu_shift;
        }
        
        $result || return;
        
        my %result = (lane => $result);
        return \%result;
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