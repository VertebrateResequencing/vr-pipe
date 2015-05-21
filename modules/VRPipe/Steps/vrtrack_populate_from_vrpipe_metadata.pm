
=head1 NAME

VRPipe::Steps::vrtrack_populate_from_vrpipe_metadata - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Steps::vrtrack_populate_from_vrpipe_metadata extends VRPipe::Steps::vrtrack_update {
    use VRPipe::Persistent::InMemory;
    use Path::Class;
    
    around options_definition {
        return {
            %{ $self->$orig },
            vrlane_storage_dir => VRPipe::StepOption->create(description => 'the absolute path to the base directory that lane-related files will be stored in (unnecessary for sequencing data)', optional => 1)
        };
    }
    
    method inputs_definition {
        return {
            files => VRPipe::StepIODefinition->create(
                type            => 'any',
                description     => 'a file that needs to be tracked',
                max_files       => -1,
                check_existence => 0
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self        = shift;
            my $options     = $self->options;
            my $db          = $options->{vrtrack_db};
            my $storage_dir = $options->{vrlane_storage_dir};
            my $req         = $self->new_requirements(memory => 500, time => 1);
            foreach my $file (@{ $self->inputs->{files} }) {
                my $path = $file->path;
                my $cmd  = "use VRPipe::Steps::vrtrack_populate_from_vrpipe_metadata; VRPipe::Steps::vrtrack_populate_from_vrpipe_metadata->populate_lane_from_file(db => q[$db], file => q[$path]";
                if ($storage_dir) {
                    $cmd .= ", storage_dir => q[$storage_dir]";
                }
                $cmd .= ");";
                $self->dispatch_vrpipecode($cmd, $req);
            }
        };
    }
    
    method description {
        return "Populate or update VRTrack lanes based on input file metadata from the irods all_with_warehouse_metadata datasource.";
    }
    
    method populate_lane_from_file (ClassName|Object $self: Str :$db!, Str|File :$file!, Str|Dir :$storage_dir?) {
        $file = VRPipe::File->get(path => $file);
        my $meta     = $file->metadata;
        my $basename = $file->basename;
        
        my $lane = $basename;
        $lane =~ s/\.gz$//;
        $lane =~ s/\.[^\.]+$//;
        
        return $self->meta_to_lane(db => $db, file => $file, meta => $meta, lane => $lane, $storage_dir ? (storage_dir => $storage_dir) : ());
    }
    
    method meta_to_lane (ClassName|Object $self: Str :$db!, HashRef :$meta!, Str :$lane!, VRPipe::File :$file!, Str|Dir :$storage_dir?) {
        my $basename = $file->basename;
        my $path     = $file->path;
        $file->add_metadata({ lane => $lane }, replace_data => 1) if $file->meta_value('lane');
        
        my $type = $file->type;
        if ($type eq 'txt' || $type eq 'bin' || $type eq 'any') {
            ($type) = $basename =~ /\.([^\.]+)$/;
        }
        $file->disconnect;
        
        if ($type =~ /bam|cram/ && $meta->{total_reads} < 1000) {
            # ignore ~empty bam files
            return 1;
        }
        
        # we can't populate vrtrack without some essential metadata
        foreach my $essential (qw(study_id study_title taxon_id)) {
            die "file $path is missing essential metadata $essential\n" unless defined $meta->{$essential};
        }
        
        # we can't represent the same sample being for 2 different projects in
        # VRTrack
        if (ref($meta->{study_id}) && @{ $meta->{study_id} } > 1) {
            die "file $path belongs to more than one study (@{$meta->{study_id}}), which VRTrack can't cope with\n";
        }
        
        my %type_to_vrtrack_type = (bam => 4, cram => 6, gtc => 7, idat => 8);
        
        my $vrtrack = $self->get_vrtrack(db => $db);
        
        # check the database to see if there are control-related notes;
        # VRTrack API makes no mention of notes, so we do this with sql!
        my $sql = qq[select note_id,note from note];
        my $sth = $vrtrack->{_dbh}->prepare($sql);
        my %notes;
        if ($sth->execute()) {
            foreach (@{ $sth->fetchall_arrayref() }) {
                $notes{ $_->[0] } = $_->[1];
            }
        }
        else {
            no warnings qw(once);
            die(sprintf('Cannot check notes table: %s', $DBI::errstr));
        }
        my $add_control_note = 0;
        if (keys %notes >= 2 && $notes{1} eq 'Control' && $notes{2} eq 'Stem cell') {
            $add_control_note = 1;
        }
        
        # the VRTrack API, even when using a transaction, doesn't have proper
        # multi-process protection and we can end up creating multiple projects
        # with the same details, which then breaks subsequent calls. So we get
        # an inmemory lock on the project to ensure we only update 1 lane per
        # project at the same time
        my $im       = VRPipe::Persistent::InMemory->new;
        my $lock_key = "vrtrack_populate.$meta->{study_id}";
        $im->block_until_locked($lock_key);
        $im->maintain_lock($lock_key);
        
        my $worked = $vrtrack->transaction(
            sub {
                # get/create the lane object
                my $vrlane;
                if (VRTrack::Lane->is_name_in_database($vrtrack, $lane, $lane)) {
                    $vrlane = VRTrack::Lane->new_by_name($vrtrack, $lane);
                }
                else {
                    $vrlane = VRTrack::Lane->create($vrtrack, $lane);
                    
                    # hierarchy name by default has special chars replaced with
                    # underscores, but it must be identical to name
                    $vrlane->hierarchy_name($lane);
                    $vrlane->update;
                }
                
                # fill out lane details we might have at this point
                $vrlane->raw_reads($meta->{total_reads}) if exists $meta->{total_reads};
                $vrlane->is_paired($meta->{is_paired_read} ? 1 : 0) if exists $meta->{is_paired_read};
                unless ($type =~ /bam|cram/) {
                    my $analysis_uuid = $meta->{analysis_uuid};
                    if ($analysis_uuid) {
                        if (ref($analysis_uuid)) {
                            $analysis_uuid = $analysis_uuid->[-1];
                        }
                        $vrlane->acc($analysis_uuid);
                        
                        if ($storage_dir && defined $meta->{irods_analysis_files}) {
                            # acc gets the latest (?) analysis_uuid, but that
                            # does not necessarily correspond to where the
                            # latest analysis files come from, so instead of
                            # having the storage path being based on the
                            # analysis_uuid as we had previously, we'll have it
                            # under the same set of subdirs as they are kept in
                            # in irods
                            my @irods_files = ref($meta->{irods_analysis_files}) ? @{ $meta->{irods_analysis_files} } : ($meta->{irods_analysis_files});
                            my $irods_file = file($storage_dir, $irods_files[0]);
                            my $storage_path;
                            if ($type eq 'idat') {
                                # since there are multiple files we store the
                                # containing dir
                                $storage_path = $irods_file->dir->stringify;
                            }
                            elsif ($type eq 'gtc') {
                                # since there should be just 1 fcr file, we
                                # store that
                                $storage_path = $irods_file->stringify;
                            }
                            $vrlane->storage_path($storage_path) if $storage_path;
                        }
                    }
                }
                if (defined $meta->{manual_qc}) {
                    if ($meta->{manual_qc}) {
                        $vrlane->npg_qc_status('pass');
                    }
                    else {
                        $vrlane->npg_qc_status('fail');
                    }
                }
                $vrlane->update;
                
                # get/create a file entry
                my $vrfile = $vrlane->get_file_by_name($basename);
                unless ($vrfile) {
                    $vrfile = $vrlane->add_file($basename);
                }
                $vrfile->type($type_to_vrtrack_type{$type});
                $vrfile->md5($meta->{md5});
                $vrfile->update;
                
                # get/create the library
                my ($library_name, $library_ssid, $library_tag_sequence);
                if ($type eq 'bam' || $type eq 'cram') {
                    $library_name = $meta->{library};
                    $library_ssid = $meta->{library_id};
                }
                else {
                    if ($type eq 'gtc') {
                        $library_name = $meta->{infinium_sample};
                    }
                    elsif ($type eq 'idat') {
                        $library_name = "$meta->{sample}.$meta->{beadchip}.$meta->{beadchip_section}";
                    }
                    $library_ssid = substr($meta->{beadchip}, -4) . substr($meta->{sample_id}, -3);
                    $library_tag_sequence = $meta->{beadchip} . '_' . $meta->{beadchip_section};
                }
                die "no library name and ssid for lane $lane\n" unless $library_name && $library_ssid;
                my $vrlibrary;
                if (VRTrack::Library->is_name_in_database($vrtrack, $library_name, $library_name)) {
                    $vrlibrary = VRTrack::Library->new_by_name($vrtrack, $library_name);
                }
                else {
                    $vrlibrary = VRTrack::Library->create($vrtrack, $library_name);
                }
                $vrlibrary->ssid($library_ssid);
                $vrlibrary->library_tag_sequence($library_tag_sequence) if $library_tag_sequence;
                $vrlibrary->update;
                
                $vrlane->library_id($vrlibrary->id);
                $vrlane->update;
                
                # get/create the project
                my $vrproject = VRTrack::Project->new_by_ssid($vrtrack, $meta->{study_id});
                my $study_title = $meta->{study_title};
                $study_title =~ s/[^\x20-\x7f]+/-/g; # title is likely to have non-ascii; replace with dashes
                unless ($vrproject) {
                    $vrproject = VRTrack::Project->create($vrtrack, $study_title);
                    $vrproject->ssid($meta->{study_id});
                    $vrproject->update;
                }
                else {
                    $vrproject->name($study_title);
                    $vrproject->update;
                }
                
                # create the study
                my $study_acc = $meta->{study_accession_number} || $meta->{study_id};
                my $vrstudy = $vrproject->get_study_by_acc($study_acc);
                unless ($vrstudy) {
                    $vrstudy = $vrproject->add_study($study_acc);
                }
                
                $vrproject->study_id($vrstudy->id);
                $vrproject->update;
                
                # get/create the sample
                my $sample_name = $meta->{sample};
                my $vrsample    = $vrproject->get_sample_by_name($sample_name);
                unless ($vrsample) {
                    $vrsample = $vrproject->add_sample($sample_name);
                }
                $vrsample->hierarchy_name($meta->{sample_supplier_name} || $meta->{sample});
                $vrsample->ssid($meta->{sample_id});
                if ($add_control_note && exists $meta->{sample_control}) {
                    $vrsample->note_id($meta->{sample_control} ? 1 : 2);
                }
                $vrsample->update;
                
                $vrlibrary->sample_id($vrsample->id);
                $vrlibrary->update;
                
                # get/create the species
                my $vrspecies = VRTrack::Species->new_by_taxon_id($vrtrack, $meta->{taxon_id});
                unless ($vrspecies) {
                    $vrspecies = VRTrack::Species->create($vrtrack, $meta->{sample_common_name}, $meta->{taxon_id});
                }
                else {
                    $vrspecies->name($meta->{sample_common_name});
                    $vrspecies->update;
                }
                
                # get/create the individual
                my $individual_name = $meta->{sample_cohort} || $meta->{sample};
                my $vrindividual = $vrsample->get_individual_by_name($individual_name);
                unless ($vrindividual) {
                    $vrindividual = $vrsample->add_individual($individual_name);
                }
                unless ($vrindividual->acc) {
                    # an individual can have multiple samples taken from it, but
                    # for historic reasons we store sample_accession_number in
                    # acc; we'll just go with the first one we see
                    $vrindividual->acc($meta->{sample_accession_number}) if defined $meta->{sample_accession_number};
                }
                my $public_name = $meta->{public_name};
                if ($public_name) {
                    # since we can have multiple sample public names per
                    # individual, we try to only set this once, preferring
                    # control sample public_name
                    if ($vrindividual->alias) {
                        if ($meta->{sample_control}) {
                            $vrindividual->alias($public_name);
                        }
                    }
                    else {
                        $vrindividual->alias($public_name);
                    }
                }
                $vrindividual->species_id($vrspecies->id);
                $vrindividual->update;
                
                $vrsample->individual_id($vrindividual->id);
                $vrsample->update;
            }
        );
        
        $im->unlock($lock_key);
        
        unless ($worked) {
            $self->throw($vrtrack->{transaction_error});
        }
    }
}

1;
