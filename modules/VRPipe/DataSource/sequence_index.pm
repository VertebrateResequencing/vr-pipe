use VRPipe::Base;
    
class VRPipe::DataSource::sequence_index with VRPipe::DataSourceRole {
    use VRPipe::Parser;
    
    has '_lanes_hash' => (is => 'rw',
                          isa => 'HashRef');
    
    method _open_source {
        my $source = file($self->source)->absolute;
        my $file = VRPipe::File->get(path => $source, type => 'txt');
        $file->e || $self->throw("sequence.index file $source does not exist!");
        return VRPipe::Parser->create('sequence_index', {file => $file});
    }
    
    method lane_fastqs (Defined :$handle, Bool :$ignore_withdrawn = 1, Str|Dir :$remote_root_dir?, Str|Dir :$local_root_dir?, Bool :$require_fastqs?, Str :$platform?, Str :$analysis_group?) {
        my $lanes_hash = $self->_lanes_hash;
        
        if ($remote_root_dir && ! $local_root_dir) {
            $self->throw("when remote_root_dir is specified, local_root_dir is required");
        }
        
        if (! defined $require_fastqs) {
            $require_fastqs = $remote_root_dir ? 0 : 1;
        }
        
        unless ($lanes_hash) {
            # because there may be multiple fastq files for a lane, and the
            # records for a lane might not be sequential, it is most efficient
            # to parse the whole file in one go and store the data
            my $pr = $handle->parsed_record;
            my %mates;
            while ($handle->next_record) {
                next if $ignore_withdrawn && $pr->[20];
                $pr->[2] || $self->throw("lane is required in column 3");
                
                if ($platform) {
                    next unless $pr->[12] =~ /$platform/;
                }
                if ($analysis_group) {
                    next unless $pr->[25] =~ /$analysis_group/;
                }
                
                my $fastq = $pr->[0] || $self->throw("fastq is required in first column");
                if ($local_root_dir) {
                    $fastq = file($local_root_dir, $fastq);
                    unless ($fastq->is_absolute) {
                        $self->throw("when local_root_dir is specified, it must result in an absolute path");
                    }
                    $fastq = $fastq->stringify;
                }
                else {
                    $fastq = file($fastq);
                    unless ($fastq->is_absolute) {
                        $self->throw("when local_root_dir is not specified, absolute paths to the fastqs must be given in column 1");
                    }
                    $fastq = $fastq->stringify;
                }
                
                my $paired = 0;
                my $mate = $pr->[19];
                if ($mate) {
                    if ($local_root_dir) {
                        $mate = file($local_root_dir, $mate)->stringify;
                    }
                    else {
                        $mate = file($mate)->stringify;
                    }
                    
                    if (exists $mates{$mate}) {
                        $paired = 2;
                    }
                    else {
                        $paired = 1;
                    }
                    $mates{$fastq} = 1;
                }
                
                my $remote_path;
                if ($remote_root_dir) {
                    if ($remote_root_dir =~ /:\/\//) {
                        $remote_path = join('/', $remote_root_dir, $fastq);
                    }
                    else {
                        $remote_path = file($remote_root_dir, $fastq)->stringify;
                    }
                }
                
                my $new_metadata = { $pr->[1] ? (expected_md5 => $pr->[1]) : (),
                                     lane => $pr->[2],
                                     study => $pr->[3],
                                     study_name => $pr->[4],
                                     center_name => $pr->[5],
                                     sample_id => $pr->[8],
                                     sample => $pr->[9],
                                     population => $pr->[10],
                                     platform => $pr->[12],
                                     library => $pr->[14],
                                     insert_size => $pr->[17],
                                     withdrawn => $pr->[20],
                                     reads => $pr->[23],
                                     bases => $pr->[24],
                                     analysis_group => $pr->[25],
                                     paired => $paired,
                                     $mate ? (mate => $mate) : (),
                                     $remote_path ? (remote_path => $remote_path) : () };
                
                my $vrfile = VRPipe::File->get(path => $fastq, type => 'fq');
                my $current_metadata = $vrfile->metadata;
                if ($current_metadata && keys %$current_metadata) {
                    foreach my $meta (qw(expected_md5 reads bases)) {
                        next unless $new_metadata->{$meta};
                        if (defined $current_metadata->{$meta} && $current_metadata->{$meta} ne $new_metadata->{$meta}) {
                            $new_metadata->{fastq_file_changed} = 1;
                            last;
                        }
                    }
                    foreach my $meta (qw(lane study study_name center_name sample_id sample population platform library insert_size analysis_group)) {
                        next unless $new_metadata->{$meta};
                        if (defined $current_metadata->{$meta} && $current_metadata->{$meta} ne $new_metadata->{$meta}) {
                            $new_metadata->{hierarchical_info_changed} = 1;
                            last;
                        }
                    }
                }
                
                $vrfile->add_metadata($new_metadata,
                                      replace_data => 0);
                
                unless ($vrfile->s) {
                    $self->throw("$fastq was in sequence.index file, but not found on disc!") if $require_fastqs;
                }
                
                push(@{$lanes_hash->{$pr->[2]}}, $fastq)
            }
            $self->_lanes_hash($lanes_hash);
        }
        
        my @lanes = sort keys %$lanes_hash;
        @lanes || return;
        
        my $lane = $lanes[0];
        my $array_ref = delete $lanes_hash->{$lane};
        my %result = (paths => $array_ref,
                      lane => $lane);
        return \%result;
    }
}

1;