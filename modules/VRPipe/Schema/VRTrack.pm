
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
    use URI::Escape;
    
    my $vrpipe_schema;
    my %pluri_type_to_order = (pluripotency => 1,     novelty => 2,     pluripotency_vs_novelty => 3,       clustering => 4,       intensity => 5);
    my %pluri_type_to_size  = (pluripotency => 'big', novelty => 'big', pluripotency_vs_novelty => 'small', clustering => 'small', intensity => 'small');
    
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
                label    => 'Lane',
                unique   => [qw(unique)],                                               # to be unique but still have correct relationships, this will need to be based on file basename
                required => [qw(lane)],
                indexed  => [qw(lane run total_reads is_paired_read qcgrind_qc_status)],
                # qcgrind_qc_status is set by qcgrind and is one of (passed failed investigate gt_pending pending); not set should be treated as pending
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
        
        my $vrlib = $vrlane->closest('VRTrack', 'Library', direction => 'incoming');
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
        
        my $vrsam = $vrlib->closest('VRTrack', 'Sample', direction => 'incoming');
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
            my $vrtax = $vrsam->closest('VRTrack', 'Taxon', direction => 'incoming');
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
        my $start = $node->closest('VRTrack', 'Lane', direction => 'incoming');
        $start ||= $node->closest('VRTrack', 'Section', direction => 'incoming');
        $start || return;
        my $start_node_id = $start->node_id;
        
        my $graph = $self->graph;
        my $db    = $graph->_global_label;
        my %nodes;
        foreach my $node ($graph->_call_vrpipe_neo4j_plugin_and_parse("/get_sequencing_hierarchy/$db/$start_node_id", namespace => 'VRTrack')) {
            my $label = $node->{label};
            bless $node, "VRPipe::Schema::VRTrack::$label";
            $nodes{ lc($label) } = $node;
        }
        
        return \%nodes;
    }
    
    method node_and_hierarchy_properties ($node) {
        my $props;
        while (my ($key, $val) = each %{ $node->properties }) {
            $props->{$key} = $val;
        }
        delete $props->{path};
        delete $props->{uuid};
        delete $props->{basename};
        
        my $h = $self->get_sequencing_hierarchy($node);
        while (my ($label, $hnode) = each %{$h}) {
            while (my ($key, $val) = each %{ $hnode->properties }) {
                $props->{"vrtrack_${label}_$key"} = $val;
            }
        }
        
        return $props;
    }
    
    # returns a hashref where keys are lc(label) and values are a qc node related
    # to the input node, where qc nodes are Bam_Stats, Genotype, Verify_Bam_ID
    # and Header_Mistakes (if they exist)
    method file_qc_nodes (Object :$node?, File|Str :$path, Str :$protocol) {
        my $nodes = $self->_call_plugin_file_qc($node, $path, $protocol);
        
        my $qc_nodes = {};
        foreach my $label (qw(bam_stats genotype verify_bam_id header_mistakes)) {
            my $node = $nodes->{$label};
            if ($node) {
                $qc_nodes->{$label} = $node;
            }
        }
        
        return $qc_nodes;
    }
    
    # to make vrtrack_metadata fast, our plugin takes a file path and returns
    # all file qc nodes along with the node for the file itself and all the
    # hierarchy nodes; file_qc_nodes() will take a subset of these, while
    # vrtrack_metadata() will convert them to a simple metadata hash; neither
    # actually needs the file node details
    method _call_plugin_file_qc ($node?, $path?, $protocol?) {
        if ($node) {
            $path     = $node->protocolless_path;
            $protocol = $node->protocol;
        }
        
        $vrpipe_schema ||= VRPipe::Schema->create('VRPipe');
        my $root = uri_escape($vrpipe_schema->protocol_to_root($protocol));
        $path = uri_escape($path);
        
        my $graph = $self->graph;
        my $db    = $graph->_global_label;
        my %nodes;
        foreach my $node ($graph->_call_vrpipe_neo4j_plugin_and_parse("/vrtrack_file_qc/$db/$root/$path", namespace => 'VRTrack')) {
            my $label = $node->{label};
            next if $label eq 'FileSystemElement';
            bless $node, "VRPipe::Schema::VRTrack::$label";
            $nodes{ lc($label) } = $node;
        }
        
        return \%nodes;
    }
    
    # return both node_and_hierarchy_properties() and file_qc_nodes() metadata
    # for a file node
    method vrtrack_metadata (Object :$node?, File|Str :$path, Str :$protocol) {
        my $meta = {};
        
        my $nodes = $self->_call_plugin_file_qc($node, $path, $protocol);
        while (my ($label, $node) = each %$nodes) {
            while (my ($key, $val) = each %{ $node->properties }) {
                next if $key eq 'uuid';
                $meta->{"vrtrack_${label}_$key"} = $val;
            }
        }
        
        return $meta;
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
    
    # this is mainly for vrtrack_qc/index.html, but it lets you get all
    # VRTrack nodes with a certain label, optionally limited to those that
    # belong to your supplied groups, studies, donors or samples (supplied as
    # node ids separated by commas). Donor and Sample nodes will have the
    # extra info that get_node_by_id_with_extra_info() adds to them.
    method nodes_of_label (Str $label!, ArrayRef[Int] :$groups?, ArrayRef[Int] :$studies?, ArrayRef[Int] :$donors?, ArrayRef[Int] :$samples?) {
        my $graph = $self->graph;
        my $db    = $graph->_global_label;
        my $args  = '';
        if ($groups) {
            $args = "?groups=" . join(',', @$groups);
            if ($studies) {
                $args .= "&studies=" . join(',', @$studies);
                $args .= "&donors=" . join(',', @$donors) if $donors;
                $args .= "&samples=" . join(',', @$samples) if $samples;
            }
        }
        my @nodes = $graph->_call_vrpipe_neo4j_plugin_and_parse("/vrtrack_nodes/$db/$label$args", namespace => 'VRTrack', label => $label);
        
        # sort by epoch and convert to ymd
        my $key = $label eq 'Donor' ? 'last_sample_added_date' : ($label eq 'Sample' ? 'created_date' : undef);
        if ($key) {
            @nodes = sort { (defined $b->{properties}->{$key} ? $b->{properties}->{$key} : 1) <=> (defined $a->{properties}->{$key} ? $a->{properties}->{$key} : 1) } @nodes;
            foreach my $node (@nodes) {
                if (defined $node->{properties}->{$key}) {
                    $node->{properties}->{$key} = DateTime->from_epoch(epoch => $node->{properties}->{$key})->ymd;
                }
            }
        }
        
        return \@nodes;
    }
    
    # adds extra info to Donor (example_sample and last_sample_added_date) and
    # Sample (donor and study ids) nodes
    method get_node_by_id_with_extra_info (Str $label!, Int $id!) {
        my $graph = $self->graph;
        my $db    = $graph->_global_label;
        return $graph->_call_vrpipe_neo4j_plugin_and_parse("/get_node_with_extra_info/$db/$id", namespace => 'VRTrack', label => $label);
    }
    
    # get all qc-related stuff associated with a donor and its samples. Also
    # get info on which samples have been QC failed/selected, along with config
    # info on who is allowed to change these and what reasons are allowed
    method donor_qc (Int $donor, Str $user, ArrayRef $new) {
        # this used to be done in multiple separate methods and queries, but now
        # it's one massive query for efficiency, and so we can do "auto qc"
        # where we add some results that summarise the specific results
        my $graph = $self->graph;
        my $db    = $graph->_global_label;
        my $args  = '';
        if ($new && @$new == 3 && $new->[0] && $new->[1] && $new->[1] =~ /^(?:failed|passed|selected|pending)$/) {
            $args = "?sample=$new->[0]&status=$new->[1]";
            if ($new->[2]) {
                $args .= "&reason=$new->[2]";
            }
        }
        my $data = $graph->_call_vrpipe_neo4j_plugin("/donor_qc/$db/$user/$donor$args");
        
        my @results;
        
        # stuff from admin details
        my $admin_details = $data->{admin_details};
        push(@results, { type => 'admin', is_admin => int($admin_details->{is_admin}), allowed_fail_reasons => $admin_details->{qc_fail_reasons} });
        
        # stuff from donor details
        my $donor_details = $data->{donor_details};
        
        # copy number by chr plot
        push(@results, { type => 'copy_number_plot', plot => $donor_details->{copy_number_by_chromosome_plot} }) if defined $donor_details->{copy_number_by_chromosome_plot};
        
        # pluritest plots
        my @pluritest_plot_results;
        foreach my $ppkey (grep { $_ =~ /^pluritest_plot_/ } keys %$donor_details) {
            my ($type) = $ppkey =~ /^pluritest_plot_(\S+)/;
            my $path = $donor_details->{$ppkey};
            next unless -s $path;
            push(@pluritest_plot_results, { type => 'pluritest_plot', path => $path, display_size => $pluri_type_to_size{$type}, order => $pluri_type_to_order{$type} });
        }
        push(@results, sort { $a->{order} <=> $b->{order} } @pluritest_plot_results);
        
        my $sample_details = $data->{samples};
        
        # first create a mapping of all sample names to some basic props, and
        # get the control sample
        my (%name_to_basic_props, $control_sample_name);
        while (my ($node_id, $sample_props) = each %$sample_details) {
            $name_to_basic_props{ $sample_props->{name} } = [$sample_props->{public_name}, $sample_props->{control}, $sample_props->{study_ids}, $node_id, $sample_props->{qc_status}];
            
            if ($sample_props->{control} == 1) {
                if (!$control_sample_name || $sample_props->{qc_status} ne 'failed') {
                    $control_sample_name = $sample_props->{public_name} . '_' . $sample_props->{name};
                }
            }
        }
        
        my (%done, %cnv_plot_paths);
        foreach my $sample_node_id (sort { $sample_details->{$b}->{control} <=> $sample_details->{$a}->{control} || $sample_details->{$a}->{public_name} cmp $sample_details->{$b}->{public_name} || $sample_details->{$a}->{name} cmp $sample_details->{$b}->{name} } keys %{$sample_details}) {
            my $sample_props     = $sample_details->{$sample_node_id};
            my $this_name        = $sample_props->{name};
            my $this_public_name = $sample_props->{public_name};
            my $this_control     = $sample_props->{control};
            my $this_study_ids   = $sample_props->{study_ids};
            my %common_results   = (sample_name => $this_name, sample_public_name => $this_public_name);
            
            # status
            push(
                @results,
                {
                    type => 'sample_status',
                    %common_results,
                    qc_status        => $sample_props->{qc_status},
                    qc_by            => $sample_props->{qc_by},
                    qc_time          => $sample_props->{qc_time} ? (DateTime->from_epoch(epoch => $sample_props->{qc_time})->ymd) : undef,
                    qc_failed_reason => $sample_props->{qc_failed_reason}
                }
            );
            
            next if $sample_props->{qc_status} eq 'failed';
            
            # gender
            push(@results, { type => 'gender', %common_results, expected_gender => $sample_props->{expected_gender}, actual_gender => $sample_props->{actual_gender}, result_file => $sample_props->{actual_gender_result_file} });
            
            # discordance
            my @disc_results;
            foreach my $type (qw(discordance_fluidigm discordance_genotyping)) {
                my $calls_string = $sample_props->{$type};
                next unless $calls_string;
                my $calls = $graph->json_decode($calls_string);
                foreach my $call (@$calls) {
                    my ($discordance, $num_of_sites, $avg_min_depth, $other_sample) = @$call;
                    next if exists $done{$type}->{$other_sample};
                    next if $this_name eq $other_sample;
                    my $other_sample_props = $name_to_basic_props{$other_sample};
                    next if $other_sample_props->[4] eq 'failed';
                    
                    my $other_study_ids = $other_sample_props->[2];
                    my $other_public    = $other_sample_props->[0];
                    if ($other_study_ids ne $this_study_ids) {
                        next unless $other_public eq $this_public_name;
                    }
                    
                    push(@disc_results, { type => $type, discordance => $discordance, num_of_sites => $num_of_sites, avg_min_depth => $avg_min_depth, sample1_name => $this_name, sample1_public_name => $this_public_name, sample1_control => $this_control, sample1_node_id => $sample_node_id, sample2_name => $other_sample, sample2_public_name => $other_public, sample2_control => $other_sample_props->[1], sample2_study => $other_study_ids, sample2_node_id => $other_sample_props->[3], study_sort => $other_study_ids eq $this_study_ids ? 1 : 2 });
                }
                $done{$type}->{$this_name} = 1;
            }
            push(@results, sort { $a->{study_sort} <=> $b->{study_sort} || $a->{sample2_study} cmp $b->{sample2_study} || $b->{sample1_control} <=> $a->{sample1_control} || $b->{sample2_control} <=> $a->{sample2_control} || $a->{sample1_public_name} cmp $b->{sample1_public_name} || $a->{sample2_public_name} cmp $b->{sample2_public_name} || $a->{sample2_name} cmp $b->{sample2_name} } @disc_results);
            
            # cnv and loh calls
            foreach my $cpkey (grep { $_ =~ /^cnv_plot_/ } keys %$sample_props) {
                my ($chr) = $cpkey =~ /^cnv_plot_(\S+)/;
                $cnv_plot_paths{$chr}->{$this_name} = $sample_props->{$cpkey};
            }
            
            my $combined_name = $this_public_name . '_' . $this_name;
            foreach my $type (qw(cnv_calls loh_calls)) {
                my $calls_string = $sample_props->{$type};
                next unless $calls_string;
                my $calls = $graph->json_decode($calls_string);
                
                if ($type eq 'cnv_calls') {
                    my $rgs = delete $calls->{RG};
                    push(@results, { type => 'copy_number_summary', sample => $combined_name, %{$calls} });
                    
                    if ($rgs) {
                        foreach my $rg (@$rgs) {
                            my $chr  = $rg->{chr};
                            my $plot = $sample_props->{"cnv_plot_$chr"};
                            push(@results, { type => 'aberrant_regions', sample => $combined_name, %{$rg}, graph => $plot ? $plot : undef });
                            $done{cnv_chrs}->{$chr} = 1;
                        }
                    }
                }
                else {
                    my @calls;
                    foreach my $call (@$calls) {
                        delete $call->{type}; # always SNPs?
                        push(@calls, { type => 'loh_calls', sample => $combined_name, control_sample => $control_sample_name, %$call });
                    }
                    push(@results, sort { $a->{sample} cmp $b->{sample} || ncmp($a->{chr}, $b->{chr}) } @calls);
                }
            }
            
            # pluritest calls
            my $calls_string = $sample_props->{pluritest_summary};
            if ($calls_string) {
                my %plur_data = %{ $graph->json_decode($calls_string) };
                my %data;
                while (my ($key, $val) = each %plur_data) {
                    $key =~ s/[- ]+/_/g;
                    $data{$key} = $val;
                }
                push(@results, { type => 'pluritest_summary', sample => $combined_name, %data });
            }
        }
        
        foreach my $chr (keys %{ $done{cnv_chrs} }) {
            delete $cnv_plot_paths{$chr};
        }
        foreach my $chr (nsort(keys %cnv_plot_paths)) {
            my $s_hash = $cnv_plot_paths{$chr};
            foreach my $sample (keys %name_to_basic_props) {
                next if $name_to_basic_props{$sample}->[4] eq 'failed';
                push(@results, { type => 'aberrant_polysomy', chr => $chr, graph => $s_hash->{$sample}, sample => $name_to_basic_props{$sample}->[0] . '_' . $sample }) if $s_hash->{$sample};
            }
        }
        
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
        my $all_samples         = $identifier eq 'sample';
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
            unless ($all_samples) {
                next if $props->{qc_failed};
            }
            $sample_meta{$sid} = [$props->{name}, $props->{public_name}, $props->{control}, $study_id, $sid, $props->{donor_node_id}, $props->{qc_failed}];
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
            $disc->{properties}->{sample} = $props->{name} unless $all_samples;
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
                if ($all_samples ? 1 : exists $allowed_samples{ $val->[3] }) {
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
            
            if ($all_samples) {
                my $sample_meta = $samples_meta[0];
                push(@$results, { type => $type, discordance => $discordance, num_of_sites => $num_of_sites, avg_min_depth => $avg_min_depth, sample_name => $sample_meta->[0], sample_public_name => $sample_meta->[1], sample_control => $sample_meta->[2] || 0, sample_node_id => $sample_meta->[4], donor_node_id => $sample_meta->[5], failed => $sample_meta->[6] || 0 });
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
