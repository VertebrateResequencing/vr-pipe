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
    
    method lane_fastqs (Defined :$handle, Bool :$ignore_withdrawn = 1, Str|Dir :$fastq_root_dir?, Bool :$require_fastqs = 1, Str :$platform?) {
        my $lanes_hash = $self->_lanes_hash;
        
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
                    next unless $pr->[12] eq $platform;
                }
                
                my $fastq = $pr->[0] || $self->throw("fastq is required in first column");
                if ($fastq_root_dir) {
                    $fastq = file($fastq_root_dir, $fastq)->absolute->stringify;
                }
                else {
                    $fastq = file($fastq)->absolute->stringify;
                }
                
                my $paired = 0;
                my $mate = $pr->[19];
                if ($mate) {
                    if ($fastq_root_dir) {
                        $mate = file($fastq_root_dir, $mate)->absolute->stringify;
                    }
                    else {
                        $mate = file($mate)->absolute->stringify;
                    }
                    
                    if (exists $mates{$mate}) {
                        $paired = 2;
                    }
                    else {
                        $paired = 1;
                    }
                    $mates{$fastq} = 1;
                }
                
                my $vrfile = VRPipe::File->get(path => $fastq, type => 'fq');
                $vrfile->add_metadata({ $pr->[1] ? (expected_md5 => $pr->[1]) : (),
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
                                        $mate ? (mate => $mate) : () },
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