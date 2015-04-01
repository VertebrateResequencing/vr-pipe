
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

class VRPipe::Schema::VRTrack with VRPipe::SchemaRole {
    use DateTime;
    use Sort::Naturally;
    my $vrpipe_schema;
    
    method schema_definitions {
        return [
            # general
            {
                label  => 'Group',  # equivalent of old mysql database name, for grouping studies that we will analyse the same way
                unique => [qw(name)]
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
                indexed      => [qw(id public_name supplier_name accession created_date consent control qc_withdrawn)],
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
            
            # QC bam stats
            {
                label          => 'Bam_Stats',
                unique         => [qw(uuid)],
                required       => [qw(mode options date), 'raw total sequences'],
                allow_anything => 1
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
    
    method add_file (Str $path) {
        $vrpipe_schema ||= VRPipe::Schema->create('VRPipe');
        return $vrpipe_schema->path_to_filesystemelement($path);
    }
    
    method get_file (Str $path) {
        $vrpipe_schema ||= VRPipe::Schema->create('VRPipe');
        return $vrpipe_schema->path_to_filesystemelement($path, only_get => 1);
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
                my $s_props = $nodes{Sample}->{$sample_id}->{properties};
                
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
            
            if (@controls || $shortest) {
                $donor->{properties}->{example_sample} = @controls == 1 ? $controls[0] : $shortest;
            }
            if ($most_recent_date) {
                # convert to yyyy-mm-dd
                my $dt = DateTime->from_epoch(epoch => $most_recent_date);
                $donor->{properties}->{last_sample_added_date} = $dt->ymd;
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
    
    # get the discordance results between all samples from the given donor.
    # these used to be stored as nodes in the graph db, but writing them all to
    # db was killing neo4j, so now we just get the gtypex file and parse it for
    # the results we need
    method donor_discordance_results (Int $donor) {
        my $study_labels = $self->cypher_labels('Study');
        
        # first we get discordance results from sequenom/fluidigm; this path is
        # specific to use of the sequenom_import_from_irods_and_covert_to_vcf +
        # vcf_merge_and_compare_genotypes pipelines
        my $cypher = "MATCH (donor)-[:sample]->(sample)-[:processed]->(irodscsv)-[:imported]->(csv)-[:converted]->(sample_vcf)-[:merged]->(merged_vcf)-[:genotypes_compared]->(gtypex) WHERE id(donor) = {donor}.id WITH gtypex MATCH (stepstate)-[sr:result]->(gtypex) WITH stepstate, sr, gtypex ORDER BY toInt(stepstate.sql_id) DESC LIMIT 1 WITH gtypex MATCH (donor)-[:sample]->(sample) WHERE id(donor) = {donor}.id WITH gtypex, sample MATCH (study:$study_labels)-[ssr:member]->(sample) RETURN DISTINCT gtypex, sample, study, ssr";
        
        my @results;
        $self->_handle_discordance_cypher($cypher, 'donor', $donor, \@results, 'discordance_fluidigm');
        
        # now get discordance results from genotyping. This path is specific to
        # use of the genome_studio_split_and_convert_to_vcf +
        # vcf_merge_and_compare_genotypes pipelines
        $cypher = "MATCH (donor)-[:sample]->(sample)-[:placed]->(section)-[:processed]->(irodsgtc)-[:imported]->(locatlgtc)-[:instigated]->(fcr)-[:converted]->(sample_vcf)-[:merged]->(merged_vcf)-[:genotypes_compared]->(gtypex) WHERE id(donor) = {donor}.id WITH gtypex MATCH (stepstate)-[sr:result]->(gtypex) WITH stepstate, sr, gtypex ORDER BY toInt(stepstate.sql_id) DESC LIMIT 1 WITH gtypex MATCH (donor)-[:sample]->(sample) WHERE id(donor) = {donor}.id WITH gtypex, sample MATCH (study:$study_labels)-[ssr:member]->(sample) RETURN DISTINCT gtypex, sample, study, ssr";
        $self->_handle_discordance_cypher($cypher, 'donor', $donor, \@results, 'discordance_genotyping');
        
        @results = sort { $a->{sample1_study} cmp $b->{sample1_study} || $a->{sample2_study} cmp $b->{sample2_study} || $b->{sample1_control} <=> $a->{sample1_control} || $b->{sample2_control} <=> $a->{sample2_control} || $a->{sample1_public_name} cmp $b->{sample1_public_name} || $a->{sample2_public_name} cmp $b->{sample2_public_name} } @results;
        
        return \@results;
    }
    
    method _handle_discordance_cypher (Str $cypher, Str $identifier, Int $node_id, ArrayRef $results, Str $type) {
        my $check_largest_study = $type eq 'discordance_fluidigm';
        my $require_sample      = $identifier eq 'sample';
        my $graph               = $self->graph;
        my $graph_data          = $graph->_run_cypher([[$cypher, { $identifier => { id => $node_id } }]]);
        
        my %rels;
        foreach my $rel (@{ $graph_data->{relationships} }) {
            push(@{ $rels{ $rel->{endNode} } }, $rel->{startNode});
        }
        
        my %nodes;
        foreach my $node (@{ $graph_data->{nodes} }) {
            $nodes{ $node->{label} }->{ $node->{id} } = $node;
        }
        
        unless (defined $nodes{Sample}) {
            return;
        }
        
        my (%sample_meta, %studies, @sample_regex, %allowed_samples);
        while (my ($sid, $sample) = each %{ $nodes{Sample} }) {
            my @study_ids;
            foreach my $study_node_id (@{ $rels{$sid} }) {
                my $study_node = $nodes{Study}->{$study_node_id};
                my $study_id   = $study_node->{properties}->{id};
                push(@study_ids, $study_id);
            }
            my $study_id = join(',', sort { $a <=> $b } @study_ids) || 0;
            $studies{$study_id}++;
            my $props = $sample->{properties};
            if ($require_sample) {
                if ($sid == $node_id) {
                    push(@sample_regex, $props->{name}, $props->{public_name});
                }
                else {
                    next if $props->{qc_withdrawn};
                }
            }
            else {
                next if $props->{qc_withdrawn};
                push(@sample_regex, $props->{name}, $props->{public_name});
            }
            $sample_meta{$sid} = [$props->{name}, $props->{public_name}, $props->{control}, $study_id, $sid];
            $allowed_samples{'public_name+sample'}->{ $props->{public_name} . '_' . $props->{name} } = $sid;
            $allowed_samples{sample}->{ $props->{name} } = $sid;
        }
        my $sample_regex = join('|', @sample_regex);
        $sample_regex = qr/$sample_regex/;
        
        my $largest_study;
        if (keys %studies > 1) {
            ($largest_study) = sort { $studies{$b} <=> $studies{$a} || $a cmp $b } keys %studies;
        }
        
        my $gtypex_path;
        while (my ($id, $node) = each %{ $nodes{FileSystemElement} }) {
            $gtypex_path = $graph->node_property($node, 'path');
            last;
        }
        return unless -s $gtypex_path;
        open(my $ifh, '<', $gtypex_path) || return;
        
        my $sample_source;
        CN: while (<$ifh>) {
            next unless /^CN\s/;
            next unless $sample_regex;
            chomp;
            my (undef, $discordance, $num_of_sites, $avg_min_depth, $sample_i, $sample_j) = split;
            
            # we're only interested in lines where both samples are samples
            # that belong to our donor and neither are qc_withdrawn.
            # complication is that the sample names in the file could be
            # name, public_name or some combination of both (or indeed
            # anything else).
            # for now we try out the obvious possibilities until found
            unless (defined $sample_source) {
                my $sample;
                if ($sample_i =~ $sample_regex) {
                    $sample = $sample_i;
                }
                else {
                    $sample = $sample_j;
                }
                
                my $sample_node;
                if ($sample =~ /^(.+)_([^_]+)$/) {
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
                    $self->throw("Couldn't find a Sample node for $sample_i in the graph database");
                }
            }
            
            my @samples_meta;
            if (defined $sample_source) {
                # there could be multiple rows with the same pair of samples; in
                # the file they are uniqufied by adding some number prefix to
                # one of the samples, but we need to get rid of that to match
                # up the names
                $sample_i =~ s/^\d+\://;
                $sample_j =~ s/^\d+\://;
                
                if ($require_sample) {
                    my $sid1 = $allowed_samples{$sample_source}->{$sample_i};
                    my $sid2 = $allowed_samples{$sample_source}->{$sample_j};
                    next CN unless ($sid1 == $node_id || $sid2 == $node_id);
                    foreach my $sid ($sid1, $sid2) {
                        push(@samples_meta, $sample_meta{$sid});
                    }
                }
                else {
                    foreach my $sample ($sample_i, $sample_j) {
                        next CN unless exists $allowed_samples{$sample_source}->{$sample};
                        push(@samples_meta, $sample_meta{ $allowed_samples{$sample_source}->{$sample} });
                    }
                }
            }
            
            if ($require_sample) {
                my $sample_meta;
                foreach my $meta (@samples_meta) {
                    next if $meta->[4] == $node_id;
                    $sample_meta = $meta;
                }
                
                push(@$results, { type => $type, discordance => $discordance, num_of_sites => $num_of_sites, avg_min_depth => $avg_min_depth, sample_name => $sample_meta->[0], sample_public_name => $sample_meta->[1], sample_control => $sample_meta->[2], sample_node_id => $sample_meta->[4] });
            }
            else {
                my ($sample1_meta, $sample2_meta) = sort { $a->[3] cmp $b->[3] || $b->[2] <=> $a->[2] || $a->[1] cmp $b->[1] } @samples_meta;
                
                if ($check_largest_study && $largest_study) {
                    # we're only interested in samples of the largest study vs
                    # all other samples in the largest study, and comparisons
                    # between samples in different studies when one of them is from
                    # the largest study and the other has the same public_name
                    if ($sample1_meta->[3] ne $largest_study || $sample2_meta->[3] ne $largest_study) {
                        if ($sample1_meta->[3] ne $largest_study && $sample2_meta->[3] ne $largest_study) {
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
        close($ifh);
    }
    
    # get the gender info for all samples from this donor
    method donor_gender_results (Int $donor) {
        my $cypher = "MATCH (donor)-[:sample]->(sample)-[ser:gender]->(egender) WHERE id(donor) = {donor}.id MATCH (sample)-[sar1:processed]->()-[sar2:imported]->()-[sar3:converted]->(resultfile)-[sar4:gender]->(agender) return sample, ser, egender, sar1, sar2, sar3, resultfile, sar4, agender";
        my $graph_data = $self->graph->_run_cypher([[$cypher, { donor => { id => $donor } }]]);
        
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
                        $result_file_path = $node->{properties}->{path};
                    }
                    
                    $previous_id = $child_node_id;
                }
                $node || next;
                
                my $gender_props = $node->{properties};
                if ($gender_props->{source} eq 'sequencescape') {
                    $eg = $gender_props->{gender};
                }
                else {
                    $ag = $gender_props->{gender};
                }
            }
            
            my $sample_props = $sample->{properties};
            push(@results, { type => 'gender', sample_name => $sample_props->{name}, sample_public_name => $sample_props->{public_name}, expected_gender => $eg, actual_gender => $ag, result_file => $result_file_path });
        }
        
        return \@results;
    }
    
    # get the discordance results between this sample and all others, along with
    # any donor info to apply to those other samples.
    # these used to be stored as nodes in the graph db, but writing them all to
    # db was killing neo4j, so now we just get the gtypex file and grep it for
    # the results we need
    method sample_discordance_results (Int $sample) {
        my $study_labels = $self->cypher_labels('Study');
        
        # first we get discordance results from sequenom/fluidigm; this path is
        # specific to use of the sequenom_import_from_irods_and_covert_to_vcf +
        # vcf_merge_and_compare_genotypes pipelines
        my $cypher = "MATCH (sample)-[:processed]->(irodscsv)-[:imported]->(csv)-[:converted]->(sample_vcf)-[:merged]->(merged_vcf)-[:genotypes_compared]->(gtypex) WHERE id(sample) = {sample}.id WITH gtypex MATCH (stepstate)-[sr:result]->(gtypex) WITH stepstate, sr, gtypex ORDER BY toInt(stepstate.sql_id) DESC LIMIT 1 WITH gtypex MATCH (allsamples)-[:processed]->(irodscsv)-[:imported]->(csv)-[:converted]->(sample_vcf)-[:merged]->(merged_vcf)-[:genotypes_compared]->(gtypex) RETURN DISTINCT allsamples, gtypex";
        
        my @results;
        $self->_handle_discordance_cypher($cypher, 'sample', $sample, \@results, 'sample_discordance_fluidigm');
        
        @results = sort { ($b->{discordance} < 3 && $b->{num_of_sites} > 15 ? 1 : 0) <=> ($a->{discordance} < 3 && $a->{num_of_sites} > 15 ? 1 : 0) || $b->{num_of_sites} - $b->{discordance} <=> $a->{num_of_sites} - $a->{discordance} || $b->{sample_control} <=> $a->{sample_control} || $a->{sample_public_name} cmp $b->{sample_public_name} } @results;
        
        return \@results;
    }
    
    method donor_cnv_results (Int $donor) {
        # the combine_bcftools_cnvs step generates a text file but doesn't
        # parse it and store the results on nodes attached to the sample,
        # because we wanted the flexibility of the file format changing and
        # being able to display anything in a table.
        # Instead we just parse the file just-in-time here when someone wants
        # the results.
        # Also get the cnv plots and copy number plot, and loh text file
        # results. We make sure to get only the latest results by getting those
        # attached to the most recently created StepState
        my $cypher = "MATCH (donor)-[:sample]->()-[:placed]->()-[:processed]->()-[:imported]->()-[:instigated]->()-[:converted]->()-[:merged]->(vcf)-[:cnv_summary|polysomy_dist]->() WHERE id(donor) = {donor}.id WITH vcf MATCH (stepstate)-[sr:result]->(vcf) WITH stepstate, sr, vcf ORDER BY toInt(stepstate.sql_id) DESC LIMIT 1 WITH vcf OPTIONAL MATCH (vcf)-[:cnv_summary]->()-[:combined]->(combined_cnvs) OPTIONAL MATCH (vcf)-[:cnv_plot]->(cnv_plots) OPTIONAL MATCH (vcf)-[:polysomy_dist]->()-[:copy_number_plot]->(copy_number_plot) RETURN DISTINCT combined_cnvs, cnv_plots, copy_number_plot ";
        # loh text file is attached to a different merged vcf
        my $cypher2 = "MATCH (donor)-[:sample]->()-[:placed]->()-[:processed]->()-[:imported]->()-[:instigated]->()-[:converted]->()-[:merged]->(vcf)-[:loh_calls]->() WHERE id(donor) = {donor}.id WITH vcf MATCH (stepstate)-[sr:result]->(vcf) WITH stepstate, sr, vcf ORDER BY toInt(stepstate.sql_id) DESC LIMIT 1 WITH vcf OPTIONAL MATCH (vcf)-[:loh_calls]->(loh_calls) RETURN DISTINCT loh_calls";
        
        my $graph       = $self->graph;
        my $graph_data  = $graph->_run_cypher([[$cypher, { donor => { id => $donor } }]]);
        my $graph_data2 = $graph->_run_cypher([[$cypher2, { donor => { id => $donor } }]]);
        
        my %cnv_plot_paths;
        my $combined_cnvs_path;
        my $copy_number_plot_path;
        my ($loh_path, $loh_control_sample);
        foreach my $node (@{ $graph_data->{nodes} }, @{ $graph_data2->{nodes} }) {
            my $basename = $graph->node_property($node, 'basename');
            my $path     = $graph->node_property($node, 'path');
            next unless -s $path;
            if ($basename eq 'combined_cnvs.txt') {
                $combined_cnvs_path = $path;
            }
            elsif ($basename eq 'copy_numbers.png') {
                $copy_number_plot_path = $path;
            }
            elsif ($basename eq 'merged.txt') {
                $loh_path = $path;
                $loh_control_sample = $graph->node_property($node, 'control_sample');
            }
            else {
                my $chr    = $graph->node_property($node, 'chr');
                my $sample = $graph->node_property($node, 'query_sample');
                $cnv_plot_paths{$chr}->{$sample} = $path;
            }
        }
        
        # parse $combined_cnvs_path
        if ($combined_cnvs_path && open(my $ifh, '<', $combined_cnvs_path)) {
            my (@samples, %results, %done_chrs, @results);
            while (<$ifh>) {
                next if /^#/;
                chomp;
                my @cols = split;
                if ($cols[0] eq 'SM') {
                    @samples = @cols[1 .. $#cols];
                }
                elsif ($cols[0] eq 'RG') {
                    my $chr = $cols[1];
                    foreach my $i (7 .. $#cols) {
                        my $sample = $samples[$i - 6];
                        push(@results, { type => 'aberrant_regions', chr => $chr, start => $cols[2], end => $cols[3], length => $cols[4], quality => $cols[5], sample => $sample, cn => $cols[$i], graph => ($cnv_plot_paths{$chr} && $cnv_plot_paths{$chr}->{$sample}) ? $cnv_plot_paths{$chr}->{$sample} : undef });
                    }
                    $done_chrs{$chr} = 1;
                }
                else {
                    foreach my $i (1 .. $#cols) {
                        $results{ $samples[$i - 1] }->{ $cols[0] } = $cols[$i];
                    }
                }
            }
            close($ifh);
            
            foreach my $chr (keys %done_chrs) {
                delete $cnv_plot_paths{$chr};
            }
            foreach my $chr (nsort(keys %cnv_plot_paths)) {
                my $s_hash = $cnv_plot_paths{$chr};
                foreach my $sample (@samples) {
                    push(@results, { type => 'aberrant_polysomy', chr => $chr, graph => $s_hash->{$sample}, sample => $sample }) if $s_hash->{$sample};
                }
            }
            
            foreach my $sample (sort { $a cmp $b } keys %results) {
                push(@results, { type => 'copy_number_summary', sample => $sample, %{ $results{$sample} } });
            }
            
            push(@results, { type => 'copy_number_plot', plot => $copy_number_plot_path });
            
            if ($loh_path && open($ifh, '<', $loh_path)) {
                my @calls;
                while (<$ifh>) {
                    chomp;
                    my ($chr, $start, $end, $sample, $count) = split(/\t/, $_);
                    push(@calls, { type => 'loh_calls', chr => $chr, start => $start, end => $end, sample => $sample, control_sample => $loh_control_sample, count => $count });
                }
                push(@results, sort { $a->{control_sample} cmp $b->{control_sample} || $a->{sample} cmp $b->{sample} || ncmp($a->{chr}, $b->{chr}) } @calls);
            }
            
            return \@results;
        }
    }
    
    method donor_pluritest_results (Int $donor) {
        # the pluritest_plot_gene_expression step generates a text file but
        # doesn't parse it and store the results on nodes attached to the
        # sample, because we wanted the flexibility of the file format changing
        # and being able to display anything in a table. Instead we just parse
        # the file just-in-time here when someone wants the results. Also get
        # the pluritest plots. We make sure to get only the latest results
        # based on coming from the most recently created stepstate
        my $cypher     = "MATCH (donor)-[:sample]->()-[:placed]->()-[:processed]->()-[:imported]->()-[:individual_profile_merge]->()-[:reformat_for_pluritest]->(reformat) WHERE id(donor) = {donor}.id WITH reformat MATCH (stepstate)-[sr:result]->(reformat) WITH stepstate, sr, reformat ORDER BY toInt(stepstate.sql_id) DESC LIMIT 1 WITH reformat OPTIONAL MATCH (reformat)-[rp:pluritest_plot]->(plots) OPTIONAL MATCH (reformat)-[rs:pluritest_summary]->(summary) RETURN distinct plots, summary";
        my $graph      = $self->graph;
        my $graph_data = $graph->_run_cypher([[$cypher, { donor => { id => $donor } }]]);
        
        my @results;
        my $plu_path;
        my %basename_to_order = ('pluritest_image02.png' => 1,     'pluritest_image03c.png' => 2,     'pluritest_image03.png' => 3,       'pluritest_image02a.png' => 4,       'pluritest_image01.png' => 5);
        my %basename_to_size  = ('pluritest_image02.png' => 'big', 'pluritest_image03c.png' => 'big', 'pluritest_image03.png' => 'small', 'pluritest_image02a.png' => 'small', 'pluritest_image01.png' => 'small');
        foreach my $node (@{ $graph_data->{nodes} }) {
            my $basename = $graph->node_property($node, 'basename');
            my $path     = $graph->node_property($node, 'path');
            next unless -s $path;
            if ($basename eq 'pluritest.csv') {
                $plu_path = $path;
            }
            else {
                push(@results, { type => 'pluritest_plot', path => $path, display_size => $basename_to_size{$basename}, order => $basename_to_order{$basename} });
            }
        }
        @results = sort { $a->{order} <=> $b->{order} } @results;
        
        # parse $plu_path
        my @summary_results;
        if ($plu_path && open(my $ifh, '<', $plu_path)) {
            <$ifh>; # header line
            while (<$ifh>) {
                chomp;
                my ($sample, $raw, $logitp, $novelty, $nov_logitp, $rmsd) = split(/,/, $_);
                $sample =~ s/^"|"$//g;
                push(@summary_results, { type => 'pluritest_summary', sample => $sample, pluri_raw => $raw, pluri_logit_p => $logitp, novelty => $novelty, novelty_logit_p => $nov_logitp, rmsd => $rmsd });
            }
            close($ifh);
        }
        push(@results, sort { $a->{sample} cmp $b->{sample} } @summary_results);
        
        return \@results if @results;
    }
}

1;
