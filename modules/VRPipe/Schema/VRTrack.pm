
=head1 NAME

VRPipe::Schema::VRTrack - schemas for tracking biology lab information

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A replacement for our vr-codebase VRTrack mysql-based database schema, this
defines all the things we need to track for our work in Vertebrate Resequencing
at the Sanger. It's pretty-much Sanger-specific.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014, 2015 Genome Research Limited.

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

class VRPipe::Schema::VRTrack with VRPipe::SchemaRole {
    use DateTime;
    use Sort::Naturally;
    my $vrpipe_schema;
    
    method schema_definitions {
        return [
            # qc-website related
            {
                label   => 'User',
                unique  => [qw(username)],
                indexed => [qw(admin)]
            },
            
            # general
            {
                label   => 'Group',              # equivalent of old mysql database name, for grouping studies that we will analyse the same way
                unique  => [qw(name)],
                indexed => [qw(qc_fail_reasons)] # we store an array of allowable reasons in qc_fail_reasons
            },
            {
                label        => 'Study',
                unique       => [qw(id)],
                indexed      => [qw(name accession)],
                keep_history => 1
            },
            {
                label   => 'Taxon',
                unique  => [qw(id)],
                indexed => [qw(common_name)]
            },
            {
                label    => 'Donor',
                unique   => [qw(id)],
                optional => [qw(example_sample last_sample_added_date)] # psuedo fields not stored in db, but you could calculate these when retrieving a donor
            },
            {
                label        => 'Sample',
                unique       => [qw(name)],
                indexed      => [qw(id public_name supplier_name accession created_date consent control qc_failed qc_selected qc_passed aberrant_chrs)], # who failed/selected/passed and for what reason is stored on a relationship between this node and a User node
                keep_history => 1
            },
            
            # bams
            {
                label        => 'Library',
                unique       => [qw(id)],
                indexed      => [qw(name tag platform center_name)],
                keep_history => 1
            },
            {
                label        => 'Lane',
                unique       => [qw(unique)],                             # to be unique but still have correct relationships, this will need to be based on file basename
                required     => [qw(lane)],
                indexed      => [qw(lane run total_reads is_paired_read)],
                keep_history => 1
            },
            {
                label  => 'Alignment',
                unique => [qw(reference)]
            },
            {
                label   => 'EBI_Run',
                unique  => [qw(acc)],
                indexed => [qw(md5)]
            },
            {
                label   => 'EBI_Submission',
                unique  => [qw(acc)],
                indexed => [qw(sub_date)]
            },
            
            # infinium
            {
                label   => 'Beadchip',
                unique  => [qw(id)],
                indexed => [qw(design)]
            },
            {
                label    => 'Section',
                unique   => [qw(unique)], # as per Lane
                required => [qw(section)]
            },
            {
                label  => 'Infinium_Sample',
                unique => [qw(id)]
            },
            {
                label  => 'Infinium_Plate',
                unique => [qw(id)]
            },
            {
                label    => 'Well',
                unique   => [qw(unique)], # as per Lane, but prefix with plate id
                required => [qw(well)]
            },
            
            # analysis files
            {
                label  => 'Analysis',
                unique => [qw(uuid)]
            },
            {
                label    => 'Collection',
                unique   => [qw(path)],
                required => [qw(date)]
            },
            
            # QC
            {
                label          => 'Bam_Stats',
                unique         => [qw(uuid)],
                required       => [qw(mode options date)],
                allow_anything => 1
            },
            {
                label          => 'Genotype',
                unique         => [qw(uuid)],
                required       => [qw(expected_sample_name matched_sample_name date pass)],
                allow_anything => 1
            },
            {
                label          => 'Verify_Bam_ID',
                unique         => [qw(uuid)],
                required       => [qw(freemix date pass)],
                allow_anything => 1
            },
            {
                label          => 'Header_Mistakes',
                unique         => [qw(uuid)],
                required       => [qw(num_mistakes)],
                indexed        => [qw(md5_of_ref_seq_md5s)],
                allow_anything => 1
            },
            {
                label    => 'Discordance',
                unique   => [qw(md5_sample)], # {md5 of the gtypex file}_{sample name}
                required => [qw(date type)],
                optional => [qw(cns)]         # this is 'optional' just so that the massive value isn't indexed
            },
            {
                label    => 'CNVs',
                unique   => [qw(md5_sample)],
                required => [qw(date)],
                optional => [qw(data)]
            },
            {
                label    => 'LOH',
                unique   => [qw(md5_sample)],
                required => [qw(date)],
                optional => [qw(data)]
            },
            {
                label    => 'Pluritest',
                unique   => [qw(md5_sample)],
                required => [qw(date)],
                optional => [qw(data)]
            },
            
            # gender (for expected or calculated from sequenom/fluidgm)
            {
                label   => 'Gender',
                unique  => [qw(source_gender_md5)],
                indexed => [qw(source gender)]
            },
        ];
    }
    
    # without enforce => 1, if the hierarchy already exists but is different,
    # no changes are made to the database and only nodes up to the difference
    # are returned
    method ensure_sequencing_hierarchy (Str :$lane!, Str :$library!, Str :$sample!, Bool :$enforce = 0, Str :$study?, Str :$taxon?, Str :$group?) {
        # get and only then try adding Lane, because we don't want to override
        # an existing lane property
        my $vrlane = $self->get('Lane', { unique => $lane });
        $vrlane ||= $self->add('Lane', { unique => $lane, lane => $lane });
        my %return = (lane => $vrlane);
        
        my ($vrlib) = $vrlane->related(incoming => { namespace => 'VRTrack', label => 'Library' });
        if ($vrlib) {
            if ($vrlib->id ne $library) {
                if ($enforce) {
                    $vrlib = $self->add('Library', { id => $library });
                    $vrlib->relate_to($vrlane, 'sequenced', selfish => 1);
                }
                else {
                    return \%return;
                }
            }
        }
        else {
            $vrlib = $self->add('Library', { id => $library }, outgoing => { type => 'sequenced', node => $vrlane });
        }
        $return{library} = $vrlib;
        
        my ($vrsam) = $vrlib->related(incoming => { namespace => 'VRTrack', label => 'Sample' });
        if ($vrsam) {
            if ($vrsam->name ne $sample) {
                if ($enforce) {
                    $vrsam = $self->add('Sample', { name => $sample });
                    $vrsam->relate_to($vrlib, 'prepared', selfish => 1);
                }
                else {
                    return \%return;
                }
            }
        }
        else {
            $vrsam = $self->add('Sample', { name => $sample }, outgoing => { type => 'prepared', node => $vrlib });
        }
        $return{sample} = $vrsam;
        
        # because a sample can be associated with more than 1 study, and study
        # is optional, we just set it if provided and don't return existing
        # studies otherwise
        if ($study) {
            $return{study} = $self->add('Study', { id => $study }, outgoing => { type => 'member', node => $vrsam });
            
            # and since a study can be in more than 1 group, we ignore group
            # if they didn't supply the study to attach it to
            if ($group) {
                $return{group} = $self->add('Group', { name => $group }, outgoing => { type => 'has', node => $return{study} });
            }
        }
        
        if ($taxon) {
            my ($vrtax) = $vrsam->related(incoming => { namespace => 'VRTrack', label => 'Taxon' });
            if ($vrtax) {
                if ($vrtax->id ne $taxon) {
                    if ($enforce) {
                        $vrtax = $self->add('Taxon', { id => $taxon });
                        $vrtax->relate_to($vrsam, 'member', selfish => 1);
                    }
                    else {
                        return \%return;
                    }
                }
            }
            else {
                $vrtax = $self->add('Taxon', { id => $taxon }, outgoing => { type => 'member', node => $vrsam });
            }
            $return{taxon} = $vrtax;
        }
        
        return \%return;
    }
    
    # given a node, we find the closest lane node pointing to it, and return the
    # lane node and the standard hierarchy of nodes above that (library, sample,
    # study, taxon, gender, donor), as a hashref keyed on lowercased node label.
    # For microarray data you could instead get (section, beadchip) instead of
    # lane and library
    method get_sequencing_hierarchy ($node, Bool :$just_preferred_study = 0) {
        my $graph        = $self->graph;
        my $taxon_labels = $self->cypher_labels('Taxon');
        my $study_labels = $self->cypher_labels('Study');
        
        my ($start) = $node->related(incoming => { max_depth => 20, namespace => 'VRTrack', label => 'Lane' });
        my $cypher_start;
        if ($start) {
            # sequencing data
            $cypher_start = 'MATCH (start)<-[:sequenced]-(second)<-[:prepared]-(sample) WHERE id(start) = {param}.start_id';
        }
        else {
            ($start) = $node->related(incoming => { max_depth => 20, namespace => 'VRTrack', label => 'Section' });
            if ($start) {
                # microarray
                $cypher_start = 'MATCH (start)<-[:placed]-(sample), (start)<-[:section]-(second) WHERE id(start) = {param}.start_id';
            }
        }
        return unless $start;
        
        # we only want to return 1 study, and will first try the query with a
        # 'preferred' study
        my $cypher = "$cypher_start OPTIONAL MATCH (sample)-[:gender]->(gender) OPTIONAL MATCH (sample)<-[:member]-(taxon:$taxon_labels) OPTIONAL MATCH (sample)<-[:sample]-(donor) OPTIONAL MATCH (sample)<-[:member { preferred: 1 }]-(study:$study_labels) RETURN start, second, sample, gender, taxon, study, donor";
        my $graph_data = $graph->_run_cypher([[$cypher, { param => { start_id => $start->node_id } }]]);
        
        my %nodes;
        foreach my $node (@{ $graph_data->{nodes} }) {
            # if there are multiple preferred study nodes, we essentially pick
            # a random one
            my $label = $node->{label};
            bless $node, "VRPipe::Schema::VRTrack::$label";
            $nodes{ lc($label) } = $node;
        }
        return unless defined $nodes{sample};
        
        unless (defined $nodes{study} || $just_preferred_study) {
            # no preferred study, we'll get and pick a random one
            $cypher = "MATCH (sample)<-[:member]-(study:$study_labels) WHERE id(sample) = {param}.sample_id RETURN study LIMIT 1";
            $graph_data = $graph->_run_cypher([[$cypher, { param => { sample_id => $nodes{sample}->node_id } }]]);
            foreach my $node (@{ $graph_data->{nodes} }) {
                my $label = $node->{label};
                bless $node, "VRPipe::Schema::VRTrack::$label";
                $nodes{ lc($label) } = $node;
            }
        }
        
        return \%nodes;
    }
    
    method add_file (Str|File $path, Str $protocol?) {
        $vrpipe_schema ||= VRPipe::Schema->create('VRPipe');
        return $vrpipe_schema->path_to_filesystemelement("$path", $protocol ? (protocol => $protocol) : ());
    }
    
    method get_file (Str|File $path, Str $protocol?) {
        $vrpipe_schema ||= VRPipe::Schema->create('VRPipe');
        if ($protocol && $protocol eq 'file:/') {
            undef $protocol;
        }
        return $vrpipe_schema->path_to_filesystemelement("$path", $protocol ? (protocol => $protocol) : (), only_get => 1);
    }
    
    # donor's just have meaningless ids; user will want to see the
    # donor's control sample public name and the most recent created_date of its
    # member samples, so we'll need some cypher that get the donor's samples
    # and the relationship between them. You're supposed to append the return
    # value of this to some cypher that first matches your desired donor(s).
    method donor_to_sample_match_cypher (Str $donor_identifer = 'donor') {
        my $sample_labels = $self->cypher_labels('Sample');
        return "MATCH ($donor_identifer)-[ds_rel]->(d_sample:$sample_labels) RETURN $donor_identifer,ds_rel,d_sample";
    }
    
    # if user ran a cypher query based on donor_to_sample_match_cypher()
    # we have sample nodes related to our donor nodes; figure out
    # which samples belong to which donors, then add the relevant
    # sample properties to the donor node. Input is the output of run_cypher().
    method add_sample_info_to_donors (HashRef $graph_data) {
        my $graph = $self->graph;
        
        my %rels;
        foreach my $rel (@{ $graph_data->{relationships} }) {
            push(@{ $rels{ $rel->{startNode} } }, $rel->{endNode});
        }
        
        my %nodes;
        foreach my $node (@{ $graph_data->{nodes} }) {
            $nodes{ $node->{label} }->{ $node->{id} } = $node;
        }
        
        foreach my $donor (values %{ $nodes{Donor} }) {
            my ($most_recent_date, @controls, $shortest);
            foreach my $sample_id (@{ $rels{ $donor->{id} } || next }) {
                my $s_props = $graph->node_properties($nodes{Sample}->{$sample_id});
                
                my $cd = $s_props->{created_date};
                if ($cd && (!defined $most_recent_date || $cd > $most_recent_date)) {
                    $most_recent_date = $cd;
                }
                
                my $name = $s_props->{public_name} || $s_props->{name};
                if (!defined $shortest || length($name) < length($shortest)) {
                    $shortest = $name;
                }
                
                push(@controls, $name) if $s_props->{control};
            }
            
            my $props = $graph->node_properties($donor);
            if (@controls || $shortest) {
                $props->{example_sample} = @controls == 1 ? $controls[0] : $shortest;
            }
            if ($most_recent_date) {
                # convert to yyyy-mm-dd
                my $dt = DateTime->from_epoch(epoch => $most_recent_date);
                $props->{last_sample_added_date} = $dt->ymd;
            }
        }
        
        return [sort { $b->{properties}->{last_sample_added_date} cmp $a->{properties}->{last_sample_added_date} || $a->{properties}->{example_sample} cmp $b->{properties}->{example_sample} } values $nodes{Donor}];
    }
    
    # it's useful when getting samples to have the donor info, if there is one,
    # so we'll need some cypher that get the sample's donor and the relationship
    # between them. You're supposed to append the return value of this to some
    # cypher that first matches your desired sample(s).
    method sample_extra_info_match_cypher (Str $sample_identifer = 'sample') {
        my $donor_labels = $self->cypher_labels('Donor');
        my $study_labels = $self->cypher_labels('Study');
        return "OPTIONAL MATCH ($sample_identifer)<-[sd_rel]-(s_donor:$donor_labels) OPTIONAL MATCH (s_study:$study_labels)-[ssr:member]->($sample_identifer) RETURN $sample_identifer,sd_rel,s_donor,ssr,s_study";
    }
    
    # if user ran a cypher query based on sample_extra_info_match_cypher()
    # we have donor and study nodes related to our sample nodes; figure out
    # which samples belong to which donors and studies, then add the relevant
    # properties to the sample node. Input is the output of run_cypher().
    method add_extra_info_to_samples (HashRef $graph_data) {
        my %rels;
        foreach my $rel (@{ $graph_data->{relationships} }) {
            push(@{ $rels{ $rel->{endNode} } }, $rel->{startNode});
        }
        
        my %nodes;
        foreach my $node (@{ $graph_data->{nodes} }) {
            $nodes{ $node->{label} }->{ $node->{id} } = $node;
        }
        
        foreach my $sample (values $nodes{Sample}) {
            my ($donor_id, $donor_node_id, $study_id, $study_node_id);
            foreach my $node_id (@{ $rels{ $sample->{id} } || next }) {
                if (exists $nodes{Donor}->{$node_id}) {
                    $donor_id      = $nodes{Donor}->{$node_id}->{properties}->{id};
                    $donor_node_id = $node_id;
                }
                elsif (exists $nodes{Study}->{$node_id}) {
                    $study_id      = $nodes{Study}->{$node_id}->{properties}->{id};
                    $study_node_id = $node_id;
                }
                last if $study_node_id && $donor_node_id;
            }
            
            if ($donor_node_id && $donor_id) {
                $sample->{properties}->{donor_node_id} = $donor_node_id;
                $sample->{properties}->{donor_id}      = $donor_id;
            }
            
            if ($study_id && $study_node_id) {
                $sample->{properties}->{study_node_id} = $study_node_id;
                $sample->{properties}->{study_id}      = $study_id;
            }
            
            # convert created_date property from epoch to ymd
            if (defined $sample->{properties}->{created_date}) {
                my $dt = DateTime->from_epoch(epoch => $sample->{properties}->{created_date});
                $sample->{properties}->{created_date} = $dt->ymd;
            }
        }
        
        return [sort { $b->{properties}->{created_date} cmp $a->{properties}->{created_date} || $a->{properties}->{public_name} cmp $b->{properties}->{public_name} } values $nodes{Sample}];
    }
    
    method get_node_by_id_with_extra_info (Str $label!, Int $id!) {
        my $graph = $self->graph();
        my $node;
        if ($label eq 'Donor') {
            # add some sample properties to the donor
            my $cypher = "MATCH (donor) WHERE id(donor) = {param}.donor_id ";
            $cypher .= $self->donor_to_sample_match_cypher('donor');
            my $graph_data = $graph->_run_cypher([[$cypher, { param => { donor_id => $id } }]]);
            my $data = $self->add_sample_info_to_donors($graph_data);
            $node = $data->[0];
        }
        elsif ($label eq 'Sample') {
            # add some donor and study properties to the sample
            my $cypher = "MATCH (sample) WHERE id(sample) = {param}.sample_id ";
            $cypher .= $self->sample_extra_info_match_cypher('sample');
            my $graph_data = $graph->_run_cypher([[$cypher, { param => { sample_id => $id } }]]);
            my $data = $self->add_extra_info_to_samples($graph_data);
            $node = $data->[0];
        }
        else {
            $self->throw("label $label not supported by node_by_id_with_extra_info");
        }
        
        return $node;
    }
    
    # get info on which samples have been QC failed/selected, along with config
    # info on who is allowed to change these and what reasons are allowed
    method donor_sample_status (Int $donor, Str $user, ArrayRef $new) {
        my $study_labels = $self->cypher_labels('Study');
        my $group_labels = $self->cypher_labels('Group');
        my $user_labels  = $self->cypher_labels('User');
        my $cypher       = "MATCH (donor)-[:sample]->(sample) WHERE id(donor) = {donor}.id MATCH (donor)<-[:member]-(study:$study_labels)<-[:has]-(group:$group_labels) OPTIONAL MATCH (group)<-[ar:administers]-(auser:$user_labels) OPTIONAL MATCH (sample)-[fbr:failed_by]->(fuser:$user_labels) OPTIONAL MATCH (sample)-[sbr:selected_by]->(suser:$user_labels) OPTIONAL MATCH (sample)-[pbr:passed_by]->(puser:$user_labels) return donor,sample,group,ar,auser,fbr,fuser,sbr,suser,pbr,puser";
        my $graph        = $self->graph;
        my $graph_data   = $graph->_run_cypher([[$cypher, { donor => { id => $donor } }]]);
        
        my %nodes = ();
        my %sample_nodes;
        my %sample_nodes_by_name;
        my %user_nodes;
        foreach my $node (@{ $graph_data->{nodes} }) {
            $nodes{ $node->{id} } = $node;
            if ($node->{label} eq 'Sample') {
                $sample_nodes{ $node->{id} } = $node;
                bless($node, 'VRPipe::Schema::VRTrack::Sample');
                $sample_nodes_by_name{ $graph->node_property($node, 'name') } = $node;
            }
            elsif ($node->{label} eq 'User') {
                $user_nodes{ $graph->node_property($node, 'username') } = $node;
            }
        }
        
        my (%failed, %selected, %passed, @fail_reasons, %fail_reasons, $is_admin);
        foreach my $rel (@{ $graph_data->{relationships} }) {
            my $end_node = $nodes{ $rel->{endNode} } || next;
            my $end_node_props = $graph->node_properties($end_node);
            if ($rel->{type} eq 'failed_by') {
                $failed{ $rel->{startNode} } = [$end_node_props->{username}, DateTime->from_epoch(epoch => $rel->{properties}->{time})->ymd, $rel->{properties}->{reason}];
            }
            elsif ($rel->{type} eq 'selected_by') {
                $selected{ $rel->{startNode} } = [$end_node_props->{username}, DateTime->from_epoch(epoch => $rel->{properties}->{time})->ymd];
            }
            elsif ($rel->{type} eq 'passed_by') {
                $passed{ $rel->{startNode} } = [$end_node_props->{username}, DateTime->from_epoch(epoch => $rel->{properties}->{time})->ymd];
            }
            elsif ($rel->{type} eq 'administers') {
                my $user_name = $graph->node_property($nodes{ $rel->{startNode} }, 'username');
                if ($user eq $user_name) {
                    $is_admin = 1;
                    foreach my $reason (@{ $end_node_props->{qc_fail_reasons} || [] }) {
                        next if exists $fail_reasons{$reason};
                        push(@fail_reasons, $reason);
                        $fail_reasons{$reason} = 1;
                    }
                }
            }
        }
        
        # set new qc status if supplied
        if ($is_admin && @$new == 3 && $new->[0] && $new->[1] && $new->[1] =~ /^(?:failed|passed|selected|pending)$/) {
            my ($sample_name, $new_status, $reason) = @$new;
            my $user_node = $user_nodes{$user};
            my $sample    = $sample_nodes_by_name{$sample_name};
            if ($sample) {
                my $time = time();
                if ($new_status eq 'selected') {
                    $sample->qc_selected(1);
                    $sample->qc_failed(0);
                    $sample->qc_passed(0);
                    $sample->relate_to($user_node, 'selected_by', replace => 1, properties => { time => $time });
                    delete $failed{ $sample->{id} };
                    $selected{ $sample->{id} } = [$user, $time];
                }
                elsif ($new_status eq 'passed') {
                    $sample->qc_passed(1);
                    $sample->qc_selected(0);
                    $sample->qc_failed(0);
                    $sample->relate_to($user_node, 'passed_by', replace => 1, properties => { time => $time });
                    delete $failed{ $sample->{id} };
                    $passed{ $sample->{id} } = [$user, $time];
                }
                elsif ($new_status eq 'failed' && $reason) {
                    $sample->qc_failed(1);
                    $sample->qc_selected(0);
                    $sample->qc_passed(0);
                    $sample->relate_to($user_node, 'failed_by', replace => 1, properties => { time => $time, reason => $reason });
                    delete $selected{ $sample->{id} };
                    $failed{ $sample->{id} } = [$user, $time, $reason];
                }
                elsif ($new_status eq 'pending') {
                    $sample->qc_selected(0);
                    $sample->qc_passed(0);
                    $sample->qc_failed(0);
                    foreach my $user ($sample->related(outgoing => { namespace => 'VRTrack', label => 'User' })) {
                        $sample->divorce_from($user);
                    }
                    delete $failed{ $sample->{id} };
                    delete $selected{ $sample->{id} };
                    delete $passed{ $sample->{id} };
                }
            }
        }
        
        my @results;
        push(@results, { type => 'admin', is_admin => $is_admin ? 1 : 0, allowed_fail_reasons => \@fail_reasons });
        
        foreach my $sample (sort { $b->{properties}->{control} <=> $a->{properties}->{control} || $a->{properties}->{public_name} cmp $b->{properties}->{public_name} || $a->{properties}->{name} cmp $b->{properties}->{name} } values %sample_nodes) {
            my $sample_props = $graph->node_properties($sample);
            my $failed       = $failed{ $sample->{id} } || [];
            my $selected     = $selected{ $sample->{id} } || [];
            my $passed       = $passed{ $sample->{id} } || [];
            push(
                @results,
                {
                    type               => 'sample_status',
                    sample_name        => $sample_props->{name},
                    sample_public_name => $sample_props->{public_name},
                    qc_status          => $sample_props->{qc_failed} ? 'failed' : ($sample_props->{qc_selected} ? 'selected' : ($sample_props->{qc_passed} ? 'passed' : 'pending')),
                    qc_by   => $sample_props->{qc_failed} ? $failed->[0] : $selected->[0] || $passed->[0],
                    qc_time => $sample_props->{qc_failed} ? $failed->[1] : $selected->[1] || $passed->[1],
                    qc_failed_reason => $failed->[2]
                }
            );
        }
        
        return \@results;
    }
    
    # get the gender info for all samples from this donor
    method donor_gender_results (Int $donor) {
        my $cypher     = "MATCH (donor)-[:sample]->(sample)-[ser:gender]->(egender) WHERE id(donor) = {donor}.id MATCH (sample)-[sar1:processed]->()-[sar2:imported]->()-[sar3:converted]->(resultfile)-[sar4:gender]->(agender) return sample, ser, egender, sar1, sar2, sar3, resultfile, sar4, agender";
        my $graph      = $self->graph;
        my $graph_data = $graph->_run_cypher([[$cypher, { donor => { id => $donor } }]]);
        
        my %rels = ();
        foreach my $rel (@{ $graph_data->{relationships} }) {
            push(@{ $rels{ $rel->{startNode} } }, $rel->{endNode});
        }
        
        my %nodes = ();
        my %sample_nodes;
        foreach my $node (@{ $graph_data->{nodes} }) {
            $nodes{ $node->{id} } = $node;
            $sample_nodes{ $node->{id} } = $node if $node->{label} eq 'Sample';
        }
        
        my @results;
        foreach my $sample (sort { $b->{properties}->{control} <=> $a->{properties}->{control} || $a->{properties}->{public_name} cmp $b->{properties}->{public_name} || $a->{properties}->{name} cmp $b->{properties}->{name} } values %sample_nodes) {
            next if $graph->node_property($sample, 'qc_failed');
            my $sid = $sample->{id};
            my ($eg, $ag, $result_file_path);
            foreach my $node_id (@{ $rels{$sid} }) {
                my $node        = $nodes{$node_id};
                my $previous_id = $node_id;
                while ($node->{label} ne 'Gender') {
                    my ($child_node_id) = @{ $rels{$previous_id} };
                    $node = $nodes{$child_node_id};
                    $node || last;
                    
                    if ($node->{label} eq 'FileSystemElement') {
                        $result_file_path = $graph->node_property($node, 'path');
                    }
                    
                    $previous_id = $child_node_id;
                }
                $node || next;
                
                my $gender_props = $graph->node_properties($node);
                if ($gender_props->{source} eq 'sequencescape') {
                    $eg = $gender_props->{gender};
                }
                else {
                    $ag = $gender_props->{gender};
                }
            }
            
            my $sample_props = $graph->node_properties($sample);
            push(@results, { type => 'gender', sample_name => $sample_props->{name}, sample_public_name => $sample_props->{public_name}, expected_gender => $eg, actual_gender => $ag, result_file => $result_file_path });
        }
        
        return \@results;
    }
    
    # get the discordance results between all samples from the given donor
    method donor_discordance_results (Int $donor) {
        my $study_labels = $self->cypher_labels('Study');
        
        # we want the most recently created Discordance node attached to each of
        # our samples, per type (fluidigm or genotype)... I can't really get my
        # head around how to construct that query, so I just get all of them
        # and will sort it out in perl
        my $cypher = "MATCH (donor)-[:sample]->(sample)-[sdr:discordance]->(disc) WHERE id(donor) = {donor}.id MATCH (study:$study_labels)-[ssr:member]->(sample) RETURN study, ssr, sample, sdr, disc";
        my $graph_data = $self->graph->_run_cypher([[$cypher, { donor => { id => $donor } }]]);
        
        my @results;
        $self->_handle_discordance_cypher($graph_data, 'donor', $donor, \@results, 'discordance_fluidigm');
        $self->_handle_discordance_cypher($graph_data, 'donor', $donor, \@results, 'discordance_genotyping');
        
        @results = sort { $a->{sample1_study} cmp $b->{sample1_study} || $a->{sample2_study} cmp $b->{sample2_study} || $b->{sample1_control} <=> $a->{sample1_control} || $b->{sample2_control} <=> $a->{sample2_control} || $a->{sample1_public_name} cmp $b->{sample1_public_name} || $a->{sample2_public_name} cmp $b->{sample2_public_name} } @results;
        
        return \@results;
    }
    
    # get the discordance results between this sample and all others
    method sample_discordance_results (Int $sample) {
        my $group_labels  = $self->cypher_labels('Group');
        my $study_labels  = $self->cypher_labels('Study');
        my $sample_labels = $self->cypher_labels('Sample');
        my $donor_labels  = $self->cypher_labels('Donor');
        
        # we also need to get info on all samples; we hope that we can limit it
        # to samples in the studies in the groups that our sample belongs to.
        # For some reason combing these in to one query is super slow, so we
        # break it in to 2
        my $cypher     = "MATCH (group:$group_labels)-[:has]->(study:$study_labels)-[:member]->(sample:$sample_labels)-[sdr:discordance]->(disc { type: 'fluidigm' }) WHERE id(sample) = {sample}.id RETURN group, sample, sdr, disc";
        my $graph      = $self->graph;
        my $graph_data = $graph->_run_cypher([[$cypher, { sample => { id => $sample } }]]);
        
        my %group_ids;
        foreach my $node (@{ $graph_data->{nodes} }) {
            next unless $node->{label} eq 'Group';
            next if $graph->node_property($node, 'name') eq 'all_studies';
            $group_ids{ $node->{id} } = 1;
        }
        my $group_ids = join(', ', sort keys %group_ids);
        
        # we're also going to get donor info so we can get the donor_node_id for
        # each sample
        $cypher = "MATCH (group:$group_labels)-[:has]->(study:$study_labels)-[:member]->(samples:$sample_labels)<-[sdr:sample]-(donor:$donor_labels) WHERE id(group) IN [$group_ids] RETURN DISTINCT samples, sdr, donor";
        my $graph_data2 = $graph->_run_cypher([[$cypher]]);
        foreach my $node (@{ $self->add_extra_info_to_samples($graph_data2) }) {
            push(@{ $graph_data->{nodes} }, $node) unless $node->{id} == $sample;
        }
        
        my @results;
        $self->_handle_discordance_cypher($graph_data, 'sample', $sample, \@results, 'sample_discordance_fluidigm');
        
        @results = sort { ($b->{discordance} < 3 && $b->{num_of_sites} > 15 ? 1 : 0) <=> ($a->{discordance} < 3 && $a->{num_of_sites} > 15 ? 1 : 0) || $b->{num_of_sites} - $b->{discordance} <=> $a->{num_of_sites} - $a->{discordance} || $b->{sample_control} <=> $a->{sample_control} || $a->{sample_public_name} cmp $b->{sample_public_name} } @results;
        
        return \@results;
    }
    
    method _handle_discordance_cypher (HashRef $graph_data, Str $identifier, Int $node_id, ArrayRef $results, Str $type) {
        my $check_largest_study = $type eq 'discordance_fluidigm';
        my $require_sample      = $identifier eq 'sample';
        my $disc_type           = $type =~ /fluidigm/ ? 'fluidigm' : 'genotype';
        my $graph               = $self->graph;
        
        my (%in_rels, %out_rels);
        foreach my $rel (@{ $graph_data->{relationships} }) {
            push(@{ $in_rels{ $rel->{endNode} } },    $rel->{startNode});
            push(@{ $out_rels{ $rel->{startNode} } }, $rel->{endNode});
        }
        
        my %nodes;
        foreach my $node (@{ $graph_data->{nodes} }) {
            $nodes{ $node->{label} }->{ $node->{id} } = $node;
        }
        
        unless (defined $nodes{Sample}) {
            return;
        }
        
        my (%sample_meta, %studies, %allowed_samples, %discs);
        while (my ($sid, $sample) = each %{ $nodes{Sample} }) {
            my @study_ids;
            foreach my $study_node_id (@{ $in_rels{$sid} }) {
                my $study_node = $nodes{Study}->{$study_node_id};
                my $study_id = $graph->node_property($study_node, 'id');
                $studies{$study_id}++;
                push(@study_ids, $study_id);
            }
            my $study_id = join(',', sort { $a <=> $b } @study_ids) || 0;
            my $props = $graph->node_properties($sample);
            if ($require_sample) {
                if ($sid != $node_id) {
                    next if $props->{qc_failed};
                }
            }
            else {
                next if $props->{qc_failed};
            }
            $sample_meta{$sid} = [$props->{name}, $props->{public_name}, $props->{control}, $study_id, $sid, $props->{donor_node_id}];
            $allowed_samples{ $props->{name} } = $sid;
            
            my @discs;
            foreach my $disc_node_id (@{ $out_rels{$sid} }) {
                my $disc = $nodes{Discordance}->{$disc_node_id};
                if ($graph->node_property($disc, 'type') eq $disc_type) {
                    push(@discs, $disc);
                }
            }
            my ($disc) = sort { $b->{properties}->{date} <=> $a->{properties}->{date} } @discs;
            next unless $disc;
            $disc->{properties}->{sample} = $props->{name} unless $require_sample;
            $discs{ $disc->{id} } = $disc;
        }
        
        my $largest_study;
        if (keys %studies > 1) {
            ($largest_study) = sort { $studies{$b} <=> $studies{$a} || $a cmp $b } keys %studies;
        }
        
        my %cns;
        foreach my $disc (values %discs) {
            my $disc_props = $graph->node_properties($disc);
            my $cns        = $graph->json_decode($disc_props->{cns});
            my $sample     = $disc_props->{sample};
            while (my ($key, $val) = each %$cns) {
                if ($require_sample ? 1 : exists $allowed_samples{ $val->[3] }) {
                    $cns{$key} = [$val->[0], $val->[1], $val->[2], sort ($val->[3], $sample ? ($sample) : ())];
                }
            }
        }
        
        foreach my $data (values %cns) {
            my ($discordance, $num_of_sites, $avg_min_depth, $sample_i, $sample_j) = @$data;
            
            my @samples_meta;
            foreach my $sample ($sample_i, $sample_j) {
                next unless $sample;
                push(@samples_meta, $sample_meta{ $allowed_samples{$sample} });
            }
            
            if ($require_sample) {
                my $sample_meta = $samples_meta[0];
                push(@$results, { type => $type, discordance => $discordance, num_of_sites => $num_of_sites, avg_min_depth => $avg_min_depth, sample_name => $sample_meta->[0], sample_public_name => $sample_meta->[1], sample_control => $sample_meta->[2] || 0, sample_node_id => $sample_meta->[4], donor_node_id => $sample_meta->[5] });
            }
            else {
                my ($sample1_meta, $sample2_meta) = sort { $a->[3] cmp $b->[3] || $b->[2] <=> $a->[2] || $a->[1] cmp $b->[1] } @samples_meta;
                
                if ($check_largest_study && $largest_study) {
                    # we're only interested in samples of the largest study vs
                    # all other samples in the largest study, and comparisons
                    # between samples in different studies when one of them is from
                    # the largest study and the other has the same public_name
                    my %sample1_studies = map { $_ => 1 } split(/,/, $sample1_meta->[3]);
                    my %sample2_studes  = map { $_ => 1 } split(/,/, $sample2_meta->[3]);
                    if (!exists $sample1_studies{$largest_study} || !exists $sample2_studes{$largest_study}) {
                        if (!exists $sample1_studies{$largest_study} && !exists $sample2_studes{$largest_study}) {
                            next;
                        }
                        unless ($sample1_meta->[1] eq $sample2_meta->[1]) {
                            next;
                        }
                    }
                }
                
                push(@$results, { type => $type, discordance => $discordance, num_of_sites => $num_of_sites, avg_min_depth => $avg_min_depth, sample1_name => $sample1_meta->[0], sample1_public_name => $sample1_meta->[1], sample1_control => $sample1_meta->[2], sample1_study => $sample1_meta->[3], sample1_node_id => $sample1_meta->[4], sample2_name => $sample2_meta->[0], sample2_public_name => $sample2_meta->[1], sample2_control => $sample2_meta->[2], sample2_study => $sample2_meta->[3], sample2_node_id => $sample2_meta->[4] });
            }
        }
    }
    
    method donor_cnv_results (Int $donor) {
        # Get the CNV summary details. Also get the cnv plots and copy number
        # plot, and loh results.
        my $cypher     = "MATCH (donor)-[:sample]->(sample)-[ssr:cnv_calls]->(summary) WHERE id(donor) = {donor}.id WITH sample,ssr,summary OPTIONAL MATCH (sample)-[scr:cnv_plot]->(cnvplot) OPTIONAL MATCH (sample)-[scnr:copy_number_by_chromosome_plot]->(cnplot) OPTIONAL MATCH (sample)-[slr:loh_calls]->(loh) RETURN sample,ssr,summary,scr,cnvplot,scnr,cnplot,slr,loh";
        my $graph      = $self->graph;
        my $graph_data = $graph->_run_cypher([[$cypher, { donor => { id => $donor } }]]);
        
        my %rels;
        foreach my $rel (@{ $graph_data->{relationships} }) {
            push(@{ $rels{ $rel->{startNode} } }, $rel->{endNode});
        }
        
        my %nodes;
        foreach my $node (@{ $graph_data->{nodes} }) {
            $nodes{ $node->{label} }->{ $node->{id} } = $node;
        }
        
        my ($copy_number_plot_path, $loh_control_sample, %cnv_plot_paths);
        while (my ($sid, $plot) = each %{ $nodes{FileSystemElement} }) {
            my $props = $graph->node_properties($plot);
            my $path  = $props->{path};
            
            my $chr          = $props->{chr};
            my $query_sample = $props->{query_sample};
            if ($chr && $query_sample) {
                # this is a sample-specific cnv plot
                $cnv_plot_paths{$chr}->{$query_sample} = $path;
                $loh_control_sample ||= $props->{control_sample}; # this isn't stored on the LOH node, so we just grab it from here instead
            }
            else {
                # this is the copy numbers plot shared by all the samples
                $copy_number_plot_path = $path;
            }
        }
        
        my @results = ({ type => 'copy_number_plot', plot => $copy_number_plot_path });
        
        my (%done_chrs, @samples);
        foreach my $sid (sort { $nodes{Sample}->{$b}->{properties}->{control} <=> $nodes{Sample}->{$a}->{properties}->{control} || $nodes{Sample}->{$a}->{properties}->{public_name} cmp $nodes{Sample}->{$b}->{properties}->{public_name} } keys %{ $nodes{Sample} }) {
            my $sample = $nodes{Sample}->{$sid};
            my $props  = $graph->node_properties($sample);
            next if $props->{qc_failed};
            my $sample_name = $props->{public_name} . '_' . $props->{name};
            push(@samples, $sample_name);
            
            foreach my $label (qw(CNVs LOH)) {
                my @nodes;
                foreach my $node_id (@{ $rels{$sid} }) {
                    push(@nodes, $nodes{$label}->{$node_id});
                }
                my ($node) = sort { ($b->{properties}->{date} || 0) <=> ($a->{properties}->{date} || 0) } @nodes;
                next unless $node;
                
                my $data = $graph->json_decode($graph->node_property($node, 'data'));
                
                if ($label eq 'CNVs') {
                    my $rgs = delete $data->{RG};
                    push(@results, { type => 'copy_number_summary', sample => $sample_name, %{$data} });
                    
                    if ($rgs) {
                        foreach my $rg (@$rgs) {
                            my $chr = $rg->{chr};
                            push(@results, { type => 'aberrant_regions', sample => $sample_name, %{$rg}, graph => ($cnv_plot_paths{$chr} && $cnv_plot_paths{$chr}->{$sample_name}) ? $cnv_plot_paths{$chr}->{$sample_name} : undef });
                            $done_chrs{$chr} = 1;
                        }
                    }
                }
                else {
                    my @calls;
                    foreach my $call (@$data) {
                        delete $call->{type}; # always SNPs?
                        push(@calls, { type => 'loh_calls', sample => $sample_name, control_sample => $loh_control_sample, %$call });
                    }
                    push(@results, sort { $a->{control_sample} cmp $b->{control_sample} || $a->{sample} cmp $b->{sample} || ncmp($a->{chr}, $b->{chr}) } @calls);
                }
            }
        }
        
        foreach my $chr (keys %done_chrs) {
            delete $cnv_plot_paths{$chr};
        }
        foreach my $chr (nsort(keys %cnv_plot_paths)) {
            my $s_hash = $cnv_plot_paths{$chr};
            foreach my $sample (@samples) {
                push(@results, { type => 'aberrant_polysomy', chr => $chr, graph => $s_hash->{$sample}, sample => $sample }) if $s_hash->{$sample};
            }
        }
        
        return \@results if @results;
    }
    
    method donor_pluritest_results (Int $donor) {
        # Get the pluitest plots attached to the donor, and the per-sample
        # pluritest details
        my $cypher     = "MATCH (donor)-[:sample]->(sample)-[sdr:pluritest]->(details) WHERE id(donor) = {donor}.id WITH donor,sample,sdr,details MATCH (donor)-[:pluritest_plot]->(plots) RETURN plots,sample,sdr,details";
        my $graph      = $self->graph;
        my $graph_data = $graph->_run_cypher([[$cypher, { donor => { id => $donor } }]]);
        
        my %rels;
        foreach my $rel (@{ $graph_data->{relationships} }) {
            push(@{ $rels{ $rel->{startNode} } }, $rel->{endNode});
        }
        
        my %nodes;
        foreach my $node (@{ $graph_data->{nodes} }) {
            $nodes{ $node->{label} }->{ $node->{id} } = $node;
        }
        
        my @results;
        my %type_to_order = (pluripotency => 1,     novelty => 2,     pluripotency_vs_novelty => 3,       clustering => 4,       intensity => 5);
        my %type_to_size  = (pluripotency => 'big', novelty => 'big', pluripotency_vs_novelty => 'small', clustering => 'small', intensity => 'small');
        while (my ($sid, $plot) = each %{ $nodes{FileSystemElement} }) {
            my $props    = $graph->node_properties($plot);
            my $basename = $props->{basename};
            my $path     = $props->{path};
            my $type     = $props->{type};
            next unless -s $path;
            push(@results, { type => 'pluritest_plot', path => $path, display_size => $type_to_size{$type}, order => $type_to_order{$type} });
        }
        @results = sort { $a->{order} <=> $b->{order} } @results;
        
        my @summary_results;
        while (my ($sid, $sample) = each %{ $nodes{Sample} }) {
            my $props = $graph->node_properties($sample);
            next if $props->{qc_failed};
            
            my @plurs;
            foreach my $plur_node_id (@{ $rels{$sid} }) {
                my $plur = $nodes{Pluritest}->{$plur_node_id};
                push(@plurs, $plur);
            }
            my ($plur) = sort { $graph->node_properties($b)->{date} <=> $graph->node_properties($a)->{date} } @plurs;
            next unless $plur;
            
            my %plur_data = %{ $graph->json_decode($graph->node_property($plur, 'data')) };
            my %data;
            while (my ($key, $val) = each %plur_data) {
                $key =~ s/[- ]+/_/g;
                $data{$key} = $val;
            }
            
            push(@summary_results, { type => 'pluritest_summary', sample => $props->{public_name} . '_' . $props->{name}, control => $props->{control} || 0, %data });
        }
        push(@results, sort { $b->{control} <=> $a->{control} || $a->{sample} cmp $b->{sample} } @summary_results);
        
        return \@results if @results;
    }
    
    method administer_qc_web_interface (Str :$user, ArrayRef[Str] :$admins?, Str :$group?, ArrayRef[Str] :$group_admins?, ArrayRef[Str] :$group_qc_fail_reasons?) {
        my %admin_users;
        # get admins
        foreach my $user_node ($self->get('User', { admin => 1 })) {
            $admin_users{ $user_node->username } = 1;
        }
        
        # if we have no admins, the first person to log in becomes admin
        unless (keys %admin_users) {
            $self->add('User', { username => $user, admin => 1 });
            $admin_users{$user} = 1;
        }
        
        # set admins
        if ($admins && exists $admin_users{$user}) {
            %admin_users = map { $_ => 1 } @$admins;
            my %done_users;
            foreach my $user_node ($self->get('User', { admin => 1 })) {
                my $name = $user_node->username();
                if (!exists $admin_users{$name}) {
                    $user_node->admin(0);
                }
                else {
                    $done_users{$name} = 1;
                }
            }
            
            foreach my $name (@$admins) {
                next if $done_users{$name};
                $self->add('User', { username => $name, admin => 1 });
            }
        }
        
        # get and set group admins and qc_fail_reasons
        my %group_admins;
        my @admin_of_groups;
        my %group_qc_fail_reasons;
        foreach my $group_node ($self->get('Group')) {
            my $group_name = $group_node->name;
            $group_admins{$group_name} = [];
            
            my $supplied = 0;
            my %supplied_group_admins;
            if ($group && $group_admins && $group eq $group_name && exists $admin_users{$user}) {
                $group_admins{$group_name} = $group_admins;
                %supplied_group_admins = map { $_ => 1 } @$group_admins;
                $supplied = 1;
            }
            
            my %already_administers;
            my $admin_of_this_group = 0;
            foreach my $user_node ($group_node->related(incoming => { namespace => 'VRTrack', label => 'User', max_depth => 1, type => 'administers' })) {
                my $user_name = $user_node->username();
                
                if ($supplied) {
                    if (!exists $supplied_group_admins{$user_name}) {
                        $group_node->divorce_from($user_node, 'administers');
                        next;
                    }
                    $already_administers{$user_name} = 1;
                }
                else {
                    push(@{ $group_admins{$group_name} }, $user_name);
                }
                
                if ($user_name eq $user) {
                    push(@admin_of_groups, $group_name);
                    $admin_of_this_group = 1;
                }
            }
            
            if ($supplied) {
                foreach my $user_name (keys %supplied_group_admins) {
                    next if exists $already_administers{$user_name};
                    next unless $user_name;
                    $self->add('User', { username => $user_name }, outgoing => { type => 'administers', node => $group_node });
                    
                    if ($user_name eq $user) {
                        push(@admin_of_groups, $group_name);
                        $admin_of_this_group = 1;
                    }
                }
            }
            
            if ($group && $group_qc_fail_reasons && $group eq $group_name && $admin_of_this_group) {
                if (@$group_qc_fail_reasons == 1 && $group_qc_fail_reasons->[0] eq '') {
                    $group_node->remove_property('qc_fail_reasons');
                }
                else {
                    $group_node->qc_fail_reasons($group_qc_fail_reasons);
                }
            }
            $group_qc_fail_reasons{$group_name} = $group_node->qc_fail_reasons() || [];
        }
        
        # provide admin-related data if user is an admin
        my @results;
        if (exists $admin_users{$user}) {
            # return the current set of admin users and groups and who admins
            # those
            push(@results, { type => 'admin', users => join(', ', sort keys %admin_users) });
            
            foreach my $group_name (sort keys %group_admins) {
                push(@results, { type => 'group_admins', group => $group_name, users => join(', ', sort @{ $group_admins{$group_name} }) });
            }
        }
        
        # provide group-admin-related data if user is admin of a group
        foreach my $group_name (sort @admin_of_groups) {
            push(@results, { type => 'group_config', group => $group_name, qc_fail_reasons => join(', ', @{ $group_qc_fail_reasons{$group_name} }) });
        }
        
        return \@results if @results;
    }
    
    # given some sample identifier sting, which might be the sample name, or
    # public_name concatenated with name, find out which
    method sample_source (Str $sample) {
        # remove possible appendage of _CTRL
        $sample =~ s/_CTRL$//;
        
        # examples are:
        # HPSI0114pf-eipl_QC1Hip-3780 => HPSI0114pf-eipl + QC1Hip-3780
        # HPSI0114i-eipl_1_EBISC_C-7595 => HPSI0114i-eipl_1 + EBISC_C-7595
        
        my ($sample_source, $sample_node);
        if ($sample =~ /^([^_]+_\d+)_(.+)$/ || $sample =~ /^(.+?)_(.+)$/) {
            # public_name+sample
            $sample_node = $self->get('Sample', { public_name => $1, name => $2 });
            if ($sample_node) {
                $sample_source = 'public_name+sample';
            }
        }
        if (!$sample_node) {
            # sample
            $sample_node = $self->get('Sample', { name => $sample });
            if ($sample_node) {
                $sample_source = 'sample';
            }
        }
        if (!$sample_node) {
            # ... give up for now
            $self->throw("Couldn't find a Sample node for $sample in the graph database");
        }
        return $sample_source;
    }
    
    # given the same input sample identifer you supplied to sample_source(),
    # get back a properties hashref that could be used to find the sample
    method sample_props_from_string (Str $sample, Str $sample_source) {
        # remove possible appendage of _CTRL
        $sample =~ s/_CTRL$//;
        
        my $sample_props;
        if ($sample_source eq 'public_name+sample') {
            # examples are:
            # HPSI0114pf-eipl_QC1Hip-3780 => HPSI0114pf-eipl + QC1Hip-3780
            # HPSI0114i-eipl_1_EBISC_C-7595 => HPSI0114i-eipl_1 + EBISC_C-7595
            my ($public_name, $name);
            if ($sample =~ /^([^_]+_\d+)_(.+)$/) {
                ($public_name, $name) = ($1, $2);
            }
            else {
                ($public_name, $name) = $sample =~ /^(.+?)_(.+)$/;
            }
            $sample_props = { public_name => $public_name, name => $name };
        }
        elsif ($sample_source eq 'sample') {
            $sample_props = { name => $sample };
        }
        return $sample_props;
    }
}

1;
