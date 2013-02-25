
=head1 NAME

VRPipe::DataSource::sequence_index - get pipeline inputs from a sequence.index

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

sequence.index files contain metadata describing fastq files, used for the 1000
genomes project and described here: L<http://www.1000genomes.org/formats>.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::DataSource::sequence_index with VRPipe::DataSourceTextRole {
    use VRPipe::Parser;
    
    method description {
        return "Use fastq files specified in a DCC sequence.index file, associating all the metadata available";
    }
    
    method source_description {
        return "The path to a DCC sequence.index file; format details can be found here: http://www.1000genomes.org/formats#IndexFiles";
    }
    
    method method_description (Str $method) {
        if ($method eq 'lane_fastqs') {
            return "An element will comprise all the fastqs for a single lane (read group - column 3), and the fastq files will have metadata from the other columns associated with them.";
        }
        if ($method eq 'sample_fastqs') {
            return "An element will comprise all the fastqs for a single sample (read group - column 10), and the fastq files will have metadata from the other columns associated with them.";
        }
        
        return '';
    }
    
    method _open_source {
        my $file = $self->source_file;
        return VRPipe::Parser->create('sequence_index', { file => $file });
    }
    
    method lane_fastqs (Defined :$handle!, Bool :$ignore_withdrawn = 1, Str|Dir :$remote_root_dir?, Str|Dir :$local_root_dir?, Bool :$require_fastqs?, Str :$platform?, Str :$analysis_group?) {
        if ($remote_root_dir && !$local_root_dir) {
            $self->throw("when remote_root_dir is specified, local_root_dir is required");
        }
        
        if (!defined $require_fastqs) {
            $require_fastqs = $remote_root_dir ? 0 : 1;
        }
        
        # because there may be multiple fastq files for a lane, and the
        # records for a lane might not be sequential, we parse the whole
        # file before making any elements
        my $pr = $handle->parsed_record;
        my $lanes_hash;
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
            my $mate   = $pr->[19];
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
                    $remote_root_dir =~ s/\/$//;
                    $remote_path = join('/', $remote_root_dir, $pr->[0]);
                }
                else {
                    $remote_path = file($remote_root_dir, $pr->[0])->stringify;
                }
            }
            
            my $new_metadata = {
                $pr->[1] ? (expected_md5 => $pr->[1]) : (),
                lane           => $pr->[2],
                study          => $pr->[3],
                study_name     => $pr->[4],
                center_name    => uc($pr->[5]),
                sample_id      => $pr->[8],
                sample         => $pr->[9],
                population     => $pr->[10],
                platform       => $pr->[12],
                library        => $pr->[14],
                insert_size    => $pr->[17],
                withdrawn      => $pr->[20],
                reads          => $pr->[23],
                bases          => $pr->[24],
                analysis_group => $pr->[25],
                paired         => $paired,
                $mate        ? (mate        => $mate)        : (),
                $remote_path ? (remote_path => $remote_path) : ()
            };
            
            my $vrfile           = VRPipe::File->create(path => $fastq, type => 'fq');
            my $current_metadata = $vrfile->metadata;
            my $changed          = 0;
            if ($current_metadata && keys %$current_metadata) {
                foreach my $meta (qw(expected_md5 reads bases)) {
                    next unless $new_metadata->{$meta};
                    if (defined $current_metadata->{$meta} && $current_metadata->{$meta} ne $new_metadata->{$meta}) {
                        $changed = 1;
                        last;
                    }
                }
                foreach my $meta (qw(lane study study_name center_name sample_id sample population platform library insert_size analysis_group)) {
                    next unless $new_metadata->{$meta};
                    if (defined $current_metadata->{$meta} && $current_metadata->{$meta} ne $new_metadata->{$meta}) {
                        $changed = 1;
                        last;
                    }
                }
            }
            
            $vrfile->add_metadata($new_metadata, replace_data => 0);
            
            unless ($vrfile->s) {
                $self->throw("$fastq was in sequence.index file, but not found on disc!") if $require_fastqs;
            }
            
            push(@{ $lanes_hash->{ $pr->[2] }->{paths} }, $fastq);
            if ($changed) {
                push(@{ $lanes_hash->{ $pr->[2] }->{changed} }, [$vrfile, $new_metadata]);
            }
        }
        
        my @element_args = ();
        my $did          = $self->_datasource_id;
        foreach my $lane (sort keys %$lanes_hash) {
            my $hash_ref = $lanes_hash->{$lane};
            my $result_hash = { paths => $hash_ref->{paths}, lane => $lane };
            push(@element_args, { datasource => $did, result => $result_hash });
            
            if ($hash_ref->{changed}) {
                my ($element) = VRPipe::DataElement->search({ datasource => $did, result => $result_hash });
                $element || next;
                foreach my $estate ($element->element_states) {
                    $estate->pipelinesetup->log_event("sequence_index DataSource will call start_from_scratch because file metadata changed", dataelement => $estate->dataelement->id);
                    $estate->start_from_scratch;
                }
                
                #*** problems happen if we start_from_scratch some of them, but
                # then get killed before updating the metadata...
                
                # only now that we've started from scratch do we we alter the
                # metadata
                foreach my $fm (@{ $hash_ref->{changed} }) {
                    my ($vrfile, $new_metadata) = @$fm;
                    $vrfile->add_metadata($new_metadata, replace_data => 1);
                }
            }
        }
        $self->_create_elements(\@element_args);
    }
    
    method sample_fastqs (Defined :$handle!, Bool :$ignore_withdrawn = 1, Str|Dir :$remote_root_dir?, Str|Dir :$local_root_dir?, Bool :$require_fastqs?, Str :$platform?, Str :$analysis_group?) {
        if ($remote_root_dir && !$local_root_dir) {
            $self->throw("when remote_root_dir is specified, local_root_dir is required");
        }
        
        if (!defined $require_fastqs) {
            $require_fastqs = $remote_root_dir ? 0 : 1;
        }
        
        # because there may be multiple fastq files for a sample, and the
        # records for a sample might not be sequential, we parse the whole
        # file before making any elements
        my $pr = $handle->parsed_record;
        my $samples_hash;
        my %mates;
        while ($handle->next_record) {
            next if $ignore_withdrawn && $pr->[20];
            $pr->[2] || $self->throw("lane is required in column 3");
            $pr->[9] || $self->throw("sample is required in column 10");
            
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
            my $mate   = $pr->[19];
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
                    $remote_root_dir =~ s/\/$//;
                    $remote_path = join('/', $remote_root_dir, $pr->[0]);
                }
                else {
                    $remote_path = file($remote_root_dir, $pr->[0])->stringify;
                }
            }
            
            my $new_metadata = {
                $pr->[1] ? (expected_md5 => $pr->[1]) : (),
                lane           => $pr->[2],
                study          => $pr->[3],
                study_name     => $pr->[4],
                center_name    => uc($pr->[5]),
                sample_id      => $pr->[8],
                sample         => $pr->[9],
                population     => $pr->[10],
                platform       => $pr->[12],
                library        => $pr->[14],
                insert_size    => $pr->[17],
                withdrawn      => $pr->[20],
                reads          => $pr->[23],
                bases          => $pr->[24],
                analysis_group => $pr->[25],
                paired         => $paired,
                $mate        ? (mate        => $mate)        : (),
                $remote_path ? (remote_path => $remote_path) : ()
            };
            
            my $vrfile           = VRPipe::File->create(path => $fastq, type => 'fq');
            my $current_metadata = $vrfile->metadata;
            my $changed          = 0;
            if ($current_metadata && keys %$current_metadata) {
                foreach my $meta (qw(expected_md5 reads bases)) {
                    next unless $new_metadata->{$meta};
                    if (defined $current_metadata->{$meta} && $current_metadata->{$meta} ne $new_metadata->{$meta}) {
                        $changed = 1;
                        last;
                    }
                }
                foreach my $meta (qw(lane study study_name center_name sample_id sample population platform library insert_size analysis_group)) {
                    next unless $new_metadata->{$meta};
                    if (defined $current_metadata->{$meta} && $current_metadata->{$meta} ne $new_metadata->{$meta}) {
                        $changed = 1;
                        last;
                    }
                }
            }
            
            $vrfile->add_metadata($new_metadata, replace_data => 0);
            
            unless ($vrfile->s) {
                $self->throw("$fastq was in sequence.index file, but not found on disc!") if $require_fastqs;
            }
            
            push(@{ $samples_hash->{ $pr->[9] }->{paths} }, $fastq);
            if ($changed) {
                push(@{ $samples_hash->{ $pr->[9] }->{changed} }, [$vrfile, $new_metadata]);
            }
        }
        
        my @element_args = ();
        my $did          = $self->_datasource_id;
        foreach my $sample (sort keys %$samples_hash) {
            my $hash_ref = $samples_hash->{$sample};
            my $result_hash = { paths => $hash_ref->{paths}, sample => $sample };
            push(@element_args, { datasource => $did, result => $result_hash });
            
            if ($hash_ref->{changed}) {
                my ($element) = VRPipe::DataElement->search({ datasource => $did, result => $result_hash });
                $element || next;
                foreach my $estate ($element->element_states) {
                    $estate->pipelinesetup->log_event("sequence_index DataSource will call start_from_scratch because file metadata changed", dataelement => $estate->dataelement->id);
                    $estate->start_from_scratch;
                }
                
                #*** problems happen if we start_from_scratch some of them, but
                # then get killed before updating the metadata...
                
                # only now that we've started from scratch do we we alter the
                # metadata
                foreach my $fm (@{ $hash_ref->{changed} }) {
                    my ($vrfile, $new_metadata) = @$fm;
                    $vrfile->add_metadata($new_metadata, replace_data => 1);
                }
            }
        }
        $self->_create_elements(\@element_args);
    }

}

1;
