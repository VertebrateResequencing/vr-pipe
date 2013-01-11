
=head1 NAME

VRPipe::DataSource::vrtrack - get pipeline inputs from a VRTrack database

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::DataSource::vrtrack with VRPipe::DataSourceRole {
    # eval these so that test suite can pass syntax check on this module when
    # VRTrack is not installed
    eval "use VRTrack::Factory;";
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
            return "An element will comprise all the bams for a single lane, and the bam files will have all relevant available metadata associated with them. The group_by_metadata option takes a '|' separated list of metadata keys by which dataelements will be grouped. e.g. group_by_metadata => 'sample|platform|library' will group all bams with the same sample, platform and library into one dataelement. Valid keys are: project, study, species, population, individual, sample, platform and library.";
        }
        elsif ($method eq 'lane_improved_bams') {
            return "An element will comprise all the bams output by the VRPipe improvement pipeline for a single lane, and the bam files will have all relevant available metadata associated with them. The group_by_metadata option takes a '|' separated list of metadata keys by which dataelements will be grouped. e.g. group_by_metadata => 'sample|platform|library' will group all bams with the same sample, platform and library into one dataelement. Valid keys are: project, study, species, population, individual, sample, platform and library.";
        }
        elsif ($method eq 'lane_fastqs') {
            return "An element will comprise all the fastqs for a single lane, and the fastq files will have all relevant available metadata associated with them. The group_by_metadata option takes a '|' separated list of metadata keys by which dataelements will be grouped. e.g. group_by_metadata => 'sample|platform|library' will group all bams with the same sample, platform and library into one dataelement. Valid keys are: project, study, species, population, individual, sample, platform and library.";
        }
        
        return '';
    }
    
    method _open_source {
        return VRTrack::Factory->instantiate(database => $self->source, mode => 'r');
    }
    
    method _has_changed {
        my $old      = $self->_changed_marker;
        my $checksum = $self->_vrtrack_lane_file_checksum;
        $self->_changed_marker($checksum);
        $old || return 1; # on first instantiation _changed_marker is undefined, defaults to changed in this case
        return ($checksum ne $old) ? 1 : 0;
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
        push(@$lane_changes, @{ VRTrack::Lane->_all_values_by_field($vrtrack_source, 'withdrawn', 'hierarchy_name') });
        
        # we care about options if supplied
        my $options = $self->options;
        foreach my $status (qw(auto_qc_status npg_qc_status qc_status gt_status)) {
            push(@$lane_changes, @{ VRTrack::Lane->_all_values_by_field($vrtrack_source, $status, 'hierarchy_name') }) if defined $options->{$status};
        }
        
        # we care about certain types of files depending on method
        my $method = $self->method;
        if ($method eq 'lane_fastqs') {
            push(@$lane_changes, @{ VRTrack::File->_all_values_by_field($vrtrack_source, 'md5', 'hierarchy_name', '(type=0 or type=1 or type=2) and latest=true') });
        }
        elsif ($method eq 'lane_bams') {
            push(@$lane_changes, @{ VRTrack::File->_all_values_by_field($vrtrack_source, 'md5', 'hierarchy_name', 'type=4 and latest=true') });
        }
        elsif ($method eq 'lane_improved_bams') {
            push(@$lane_changes, @{ VRTrack::File->_all_values_by_field($vrtrack_source, 'md5', 'hierarchy_name', 'type=5 and latest=true') });
        }
        
        my $digest = md5_hex join('', map { defined $_ ? $_ : 'NULL' } @$lane_changes);
        return $digest;
    }
    
    method _filtered_lanes (Defined :$handle!, Str :$project_regex?, Str :$sample_regex?, Str :$library_regex?, Str :$gt_status?, Str :$qc_status?, Str :$auto_qc_status?, Str :$npg_qc_status?) {
        my @lanes = $handle->get_lanes($project_regex ? (project_regex => $project_regex) : (), $sample_regex ? (sample_regex => $sample_regex) : (), $library_regex ? (library_regex => $library_regex) : ());
        
        my @filtered;
        foreach my $lane (@lanes) {
            if (defined $gt_status) {
                my $this_status = $lane->genotype_status;
                next if $this_status !~ /$gt_status/;
            }
            if (defined $qc_status) {
                my $this_status = $lane->qc_status;
                next if $this_status !~ /$qc_status/;
            }
            if (defined $auto_qc_status) {
                my $this_status = $lane->auto_qc_status;
                next if $this_status !~ /$auto_qc_status/;
            }
            if (defined $npg_qc_status) {
                my $this_status = $lane->npg_qc_status;
                next if $this_status !~ /$npg_qc_status/;
            }
            
            push(@filtered, $lane);
        }
        
        return @filtered;
    }
    
    method lanes (Defined :$handle!, Str :$project_regex?, Str :$sample_regex?, Str :$library_regex?, Str :$gt_status?, Str :$qc_status?, Str :$auto_qc_status?, Str :$npg_qc_status?) {
        my %args;
        $args{handle}         = $handle         if defined($handle);
        $args{project_regex}  = $project_regex  if defined($project_regex);
        $args{sample_regex}   = $sample_regex   if defined($sample_regex);
        $args{library_regex}  = $library_regex  if defined($library_regex);
        $args{gt_status}      = $gt_status      if defined($gt_status);
        $args{qc_status}      = $qc_status      if defined($qc_status);
        $args{auto_qc_status} = $auto_qc_status if defined($auto_qc_status);
        $args{npg_qc_status}  = $npg_qc_status  if defined($npg_qc_status);
        
        my @element_args;
        foreach my $lane ($self->_filtered_lanes(%args)) {
            push(@element_args, { datasource => $self->_datasource_id, result => { lane => $lane->hierarchy_name } });
        }
        $self->_create_elements(\@element_args);
    }
    
    method lane_bams (Defined :$handle!, Str|Dir :$local_root_dir!, Str :$project_regex?, Str :$sample_regex?, Str :$library_regex?, Str :$gt_status?, Str :$qc_status?, Str :$auto_qc_status?, Str :$npg_qc_status?, Str :$group_by_metadata?) {
        my %args;
        $args{handle}            = $handle            if defined($handle);
        $args{local_root_dir}    = $local_root_dir    if defined($local_root_dir);
        $args{project_regex}     = $project_regex     if defined($project_regex);
        $args{sample_regex}      = $sample_regex      if defined($sample_regex);
        $args{library_regex}     = $library_regex     if defined($library_regex);
        $args{gt_status}         = $gt_status         if defined($gt_status);
        $args{qc_status}         = $qc_status         if defined($qc_status);
        $args{auto_qc_status}    = $auto_qc_status    if defined($auto_qc_status);
        $args{npg_qc_status}     = $npg_qc_status     if defined($npg_qc_status);
        $args{group_by_metadata} = $group_by_metadata if defined($group_by_metadata);
        
        # add to the argument list to filter on bam files
        $args{'file_type'} = 4;
        return $self->_lane_files(%args);
    }
    
    method lane_improved_bams (Defined :$handle!, Str :$project_regex?, Str :$sample_regex?, Str :$library_regex?, Str :$gt_status?, Str :$qc_status?, Str :$auto_qc_status?, Str :$npg_qc_status?, Str :$group_by_metadata?) {
        my %args;
        $args{handle}            = $handle            if defined($handle);
        $args{project_regex}     = $project_regex     if defined($project_regex);
        $args{sample_regex}      = $sample_regex      if defined($sample_regex);
        $args{library_regex}     = $library_regex     if defined($library_regex);
        $args{gt_status}         = $gt_status         if defined($gt_status);
        $args{qc_status}         = $qc_status         if defined($qc_status);
        $args{auto_qc_status}    = $auto_qc_status    if defined($auto_qc_status);
        $args{npg_qc_status}     = $npg_qc_status     if defined($npg_qc_status);
        $args{group_by_metadata} = $group_by_metadata if defined($group_by_metadata);
        
        # add to the argument list to filter on improved bam files
        $args{'file_type'} = 5;
        return $self->_lane_files(%args);
    }
    
    method lane_fastqs (Defined :$handle!, Str|Dir :$local_root_dir!, Str :$project_regex?, Str :$sample_regex?, Str :$library_regex?, Str :$gt_status?, Str :$qc_status?, Str :$auto_qc_status?, Str :$npg_qc_status?, Str :$group_by_metadata?) {
        my %args;
        $args{handle}            = $handle            if defined($handle);
        $args{local_root_dir}    = $local_root_dir    if defined($local_root_dir);
        $args{project_regex}     = $project_regex     if defined($project_regex);
        $args{sample_regex}      = $sample_regex      if defined($sample_regex);
        $args{library_regex}     = $library_regex     if defined($library_regex);
        $args{gt_status}         = $gt_status         if defined($gt_status);
        $args{qc_status}         = $qc_status         if defined($qc_status);
        $args{auto_qc_status}    = $auto_qc_status    if defined($auto_qc_status);
        $args{npg_qc_status}     = $npg_qc_status     if defined($npg_qc_status);
        $args{group_by_metadata} = $group_by_metadata if defined($group_by_metadata);
        
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
        my $local_root_dir  = delete $args{local_root_dir};
        unless ($file_type eq "5") {
            $self->throw("local_root_dir is required") unless $local_root_dir;
        }
        my $vrtrack = $args{handle} || $self->throw("handle is required");
        my $group_by_metadata = delete $args{group_by_metadata};
        
        my @single_results;
        foreach my $lane ($self->_filtered_lanes(%args)) {
            my %lane_info = $vrtrack->lane_info($lane);
            my @lane_changed_details;
            
            my @files;
            foreach my $file (@{ $lane->files }) {
                next unless $file->type =~ /^($file_type)$/;
                
                my ($file_abs_path, $vrfile);
                if ($file_type eq "5") {
                    my $file_name = $file->name;
                    my ($fid) = $file_name =~ /VRPipe::File::(\d+)/;
                    $fid || $self->throw("file " . $file->id . " was type 5, but did not have a VRPipe::File name (was '$file_name')");
                    $vrfile = VRPipe::File->get(id => $fid);
                    $file_abs_path = $vrfile->path->stringify;
                }
                else {
                    $file_abs_path = file($local_root_dir, $file->name)->stringify;
                    $vrfile = VRPipe::File->create(path => $file_abs_path, type => $vrpipe_filetype);
                }
                
                my $new_metadata = {
                    expected_md5 => $file->md5,
                    center_name  => $lane_info{centre},
                    project      => $lane_info{project},
                    study        => $lane_info{study},
                    $lane_info{species} ? (species => $lane_info{species}) : (),
                    population  => $lane_info{population},
                    individual  => $lane_info{individual},
                    sample      => $lane_info{sample},
                    platform    => $lane_info{seq_tech},
                    library     => $lane_info{library},
                    lane        => $lane_info{lane},
                    withdrawn   => $lane_info{withdrawn} || 0,   #*** we don't actually handle withdrawn files properly atm; if all withdrawn we shouldn't create the element...
                    insert_size => $lane_info{insert_size} || 0,
                    reads       => $file->raw_reads || 0,
                    bases       => $file->raw_bases || 0,
                    paired      => $lane_info{vrlane}->is_paired,
                    lane_id     => $file->lane_id
                };
                
                # add metadata to file but ensure that we update any fields in
                # the new metadata
                my $current_metadata = $vrfile->metadata;
                my $changed          = 0;
                if ($current_metadata && keys %$current_metadata) {
                    foreach my $meta (qw(expected_md5 lane project study center_name sample population platform library insert_size)) {
                        next unless defined $new_metadata->{$meta};
                        if (defined $current_metadata->{$meta} && $current_metadata->{$meta} ne $new_metadata->{$meta}) {
                            $self->debug("metadata '$meta' changed from $current_metadata->{$meta} to $new_metadata->{$meta} for file $file_abs_path, so will mark lane " . $lane->name . " as changed");
                            $changed = 1;
                            last;
                        }
                    }
                    
                    # some fields we'll just blindly update the metadata for,
                    # since they do not appear in bam headers; there's no need
                    # to start_from_scratch when they change??...
                    foreach my $meta (qw(species individual)) {
                        next unless defined $new_metadata->{$meta};
                        if (defined $current_metadata->{$meta} && $current_metadata->{$meta} ne $new_metadata->{$meta}) {
                            $vrfile->add_metadata({ $meta => $new_metadata->{$meta} }, replace_data => 1);
                        }
                    }
                }
                
                # if there was no metadata this will add metadata to the file.
                $vrfile->add_metadata($new_metadata, replace_data => 0);
                
                # if there was a change in VRPipe::File metadata store it in a
                # hash and change the metadata in the VRPipe::File later when
                # more appropriate, having made sure the DataElement element
                # states have been reset (see below)
                if ($changed) {
                    push(@lane_changed_details, [$vrfile, $new_metadata]);
                }
                
                push @files, $file_abs_path;
            }
            
            push(
                @single_results,
                {
                    paths               => \@files,
                    lane                => $lane->name,
                    groupable_lane_meta => {
                        project => $lane_info{project},
                        study   => $lane_info{study},
                        $lane_info{species} ? (species => $lane_info{species}) : (),
                        population => $lane_info{population},
                        individual => $lane_info{individual},
                        sample     => $lane_info{sample},
                        platform   => $lane_info{seq_tech},
                        library    => $lane_info{library}
                    },
                    scalar(@lane_changed_details) ? (changed => \@lane_changed_details) : ()
                }
            );
        }
        
        my $results;
        my $group_name;
        if ($group_by_metadata) {
            my @meta_keys = split /\|/, $group_by_metadata;
            
            my $group_hash;
            my %groups_missing_files;
            foreach my $result (@single_results) {
                my @group_keys;
                foreach my $key (@meta_keys) {
                    $self->throw("Metadata key $key not present for lane " . $result->{lane}) unless (exists $result->{groupable_lane_meta}->{$key});
                    push @group_keys, $result->{groupable_lane_meta}->{$key};
                }
                
                my $group_key = join '|', @group_keys;
                push(@{ $group_hash->{$group_key}->{paths} }, @{ $result->{paths} });
                push(@{ $group_hash->{$group_key}->{changed} }, @{ $result->{changed} }) if $result->{changed};
                unless (@{ $result->{paths} }) {
                    $groups_missing_files{$group_key} = 1;
                }
            }
            
            while (my ($group, $data) = each %$group_hash) {
                next if exists $groups_missing_files{$group};
                push(@$results, { paths => $data->{paths}, group => $group, $data->{changed} ? (changed => $data->{changed}) : () });
            }
            
            $group_name = 'group';
        }
        else {
            $results    = \@single_results;
            $group_name = 'lane';
        }
        
        my $did = $self->_datasource_id;
        my @element_args;
        foreach my $result (@$results) {
            my $result_hash = { paths => $result->{paths}, $group_name => $result->{$group_name} };
            push(@element_args, { datasource => $did, result => $result_hash });
            
            if ($result->{changed}) {
                # because 'changed' is based on file metadata changing, and the
                # file may have had its metadata applied in some other pipeline
                # for some other datasource, there may be no DataElement to
                # actually change
                my ($element) = VRPipe::DataElement->search({ datasource => $did, result => $result_hash });
                $element || next;
                
                # reset element states first
                foreach my $estate ($element->element_states) {
                    $estate->start_from_scratch;
                }
                # then change metadata in files
                foreach my $fm (@{ $result->{changed} }) {
                    my ($vrfile, $new_metadata) = @$fm;
                    $vrfile->add_metadata($new_metadata, replace_data => 1);
                }
            }
        }
        $self->_create_elements(\@element_args);
    }
}

1;
