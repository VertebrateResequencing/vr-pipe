
=head1 NAME

VRPipe::Steps::update_vrtrack_from_graph_db - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::update_vrtrack_from_graph_db extends VRPipe::Steps::vrtrack_update {
    use VRPipe::Parser;
    use VRPipe::Schema;
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                description => 'bam files',
                max_files   => -1,
                metadata    => { lane => 'lane name (a unique identifer for this sequencing run, aka read group)' }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            my $opts = $self->options;
            my $req  = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $file (@{ $self->inputs->{bam_files} }) {
                my $path = $file->path;
                my $lane = $file->metadata->{lane};
                
                unless ($lane) {
                    # modern steps don't store stuff in ->metadata(); check the
                    # graph db instead
                    my $schema = VRPipe::Schema->create("VRPipe");
                    my $file_graph_node = $schema->get('File', { path => $file->path->stringify });
                    if ($file_graph_node) {
                        my $lane_node = $file_graph_node->closest('VRTrack', 'Lane', direction => 'incoming');
                        if ($lane_node) {
                            $lane = $lane_node->unique;
                        }
                    }
                }
                unless ($lane) {
                    $self->throw("file $path lacks lane metadata");
                }
                
                my $cmd = "use VRPipe::Steps::update_vrtrack_from_graph_db; VRPipe::Steps::update_vrtrack_from_graph_db->update_vrtrack(db => q[$opts->{vrtrack_db}], bam => q[$path], lane => q[$lane] );";
                $self->dispatch_vrpipecode($cmd, $req);
            }
        };
    }
    
    method outputs_definition {
        return {};
    }
    
    method description {
        return "Update the VRTrack database using the auto qc status for a lane in the graph db.";
    }
    
    method update_vrtrack (ClassName|Object $self: Str :$db!, Str|File :$bam!, Str :$lane! ) {
        my $bam_file = VRPipe::File->get(path => $bam);
        my $meta = $bam_file->metadata;
        
        # get the auto qc results from graph db
        my $schema = VRPipe::Schema->create('VRTrack');
        
        my $file_graph_node = $schema->get_file($bam_file->protocolless_path, $bam_file->protocol);
        unless ($file_graph_node) {
            $self->throw($bam_file->path . " was not in the graph database");
        }
        my ($auto_qc_node) = $file_graph_node->related(outgoing => { type => 'auto_qc_status' });
        my $auto_qc = $auto_qc_node->properties;
        
        # get the lane object from VRTrack
        my $vrtrack = $self->get_vrtrack(db => $db);
        my $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane) || $self->throw("No lane named '$lane' in database '$db'");
        my $mapstats = $vrlane->latest_mapping || $self->throw("There were no mapstats for lane $lane");
        
        # high level qc based on insert sizes
        my $lib_status = ($auto_qc->{'Insert size'}[0] && $auto_qc->{'Insert size (rev)'}[0]) ? 'passed' : 'failed';
        my $lib_to_update = VRTrack::Library->new_by_field_value($vrtrack, 'library_id', $vrlane->library_id()) or $self->throw("No vrtrack library for lane $lane?");
        
        # genotype check results from metadata
        my @gtype_results;
        my $gstatus;
        my $gtype_analysis = $meta->{gtype_analysis};
        if ($gtype_analysis) {
            ($gstatus) = $gtype_analysis =~ /status=(\S+) expected=(\S+) found=(\S+) (?:ratio|concordance)=(\S+)/;
            @gtype_results = ($gstatus, $2, $3, $4);
        }
        
        # write results to the VRTrack database
        my @objs_to_check = ($vrlane);
        push(@objs_to_check, $lib_to_update) if $lib_to_update;
        my $worked = $vrtrack->transaction(
            sub {
                if ($lib_to_update) {
                    # don't pass the library if another lane has previously set it
                    # to failed
                    my $do_update      = 1;
                    my $current_status = $lib_to_update->auto_qc_status;
                    if ($lib_status eq 'passed' && $current_status && $current_status eq 'failed') {
                        $do_update = 0;
                    }
                    
                    if ($do_update) {
                        $lib_to_update->auto_qc_status($lib_status);
                        $lib_to_update->update();
                    }
                }
                
                # output autoQC results to the mapping stats
                for my $key (keys %$auto_qc) {
                    next if ($key eq "pass" || $key eq "date" || $key eq "uuid");
                    $mapstats->add_autoqc($key, $auto_qc->{$key}[0], $auto_qc->{$key}[1]);
                }
                $mapstats->update;
                
                $vrlane->auto_qc_status($meta->{auto_qc_status});
                
                # also, if we did our own genotype check, write those results back
                # to VRTrack now
                if (@gtype_results) {
                    my $mapstats = $vrlane->latest_mapping;
                    $mapstats->genotype_expected($gtype_results[1]);
                    $mapstats->genotype_found($gtype_results[2]);
                    $mapstats->genotype_ratio($gtype_results[3]);
                    $mapstats->update();
                    
                    $vrlane->genotype_status($gtype_results[0]);
                }
                
                $vrlane->update();
            },
            undef,
            [@objs_to_check]
        );
        
        # for some bizarre reason, at this point $lib_to_update->auto_qc_status
        # can report the desired status, yet the database has not actually been
        # updated. Check this
        if ($worked && $lib_to_update) {
            $vrtrack = $self->get_vrtrack(db => $db);
            my $lib_id            = $lib_to_update->id;
            my $check_lib         = VRTrack::Library->new($vrtrack, $lib_id);
            my $desired_qc_status = $lib_to_update->auto_qc_status;
            my $actual_qc_status  = $check_lib->auto_qc_status;
            $self->throw("the auto_qc_status we set ($desired_qc_status) does not match the one in the db ($actual_qc_status) for lane $lib_id") unless $actual_qc_status eq $desired_qc_status;
            
            # below commented section definitely solves the problem, but latest
            # VRTrack has a more generic solution (not yet confirmed effective)
            
            #my $max_retries = 10;
            #while ($check_lib->auto_qc_status ne $desired_qc_status) {
            #    warn "library auto_qc_status in the database was not $desired_qc_status, will try and set it again...\n";
            #    $vrtrack->transaction(sub {
            #        $check_lib->auto_qc_status($desired_qc_status);
            #        $check_lib->update;
            #    });
            #
            #    $max_retries--;
            #    if ($max_retries <= 0) {
            #        $self->throw("Could not get library auto_qc_status to update in the database for library $lib_id");
            #    }
            #
            #    $vrtrack = $self->get_vrtrack(db => $db);
            #    $check_lib = VRTrack::Library->new($vrtrack, $lib_id);
            #}
            #warn "Pretty sure that library auto_qc_status in the database is now $desired_qc_status\n";
        }
        
        $self->throw($vrtrack->{transaction_error}) unless ($worked);
    }

}

1;
